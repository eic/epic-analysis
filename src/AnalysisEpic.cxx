// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek, Christopher Dilks

#include "AnalysisEpic.h"

AnalysisEpic::AnalysisEpic(TString infileName_, TString outfilePrefix_)
  : Analysis(infileName_, outfilePrefix_)
{};

AnalysisEpic::~AnalysisEpic() {};

void AnalysisEpic::Execute()
{
  // setup
  Prepare();

  // read EventEvaluator tree
  auto chain = std::make_unique<TChain>("events");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      // std::cout << "Adding " << infiles[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str());
    }
  }
  chain->CanDeleteRefs();
  auto listOfBranches = chain->GetListOfBranches();

  TTreeReader tr(chain.get());

  TTreeReaderArray<Int_t> hepmcp_status(tr, "GeneratedParticles.type");
  TTreeReaderArray<Int_t> hepmcp_PDG(tr,    "GeneratedParticles.PDG");
  TTreeReaderArray<Float_t> hepmcp_E(tr,      "GeneratedParticles.energy");
  TTreeReaderArray<Float_t> hepmcp_psx(tr,    "GeneratedParticles.momentum.x");
  TTreeReaderArray<Float_t> hepmcp_psy(tr,    "GeneratedParticles.momentum.y");
  TTreeReaderArray<Float_t> hepmcp_psz(tr,    "GeneratedParticles.momentum.z");

  

  // All true particles (including secondaries, etc)
  TTreeReaderArray<Int_t> mcpart_PDG(tr,      "MCParticles.PDG");
  TTreeReaderArray<Int_t> mcpart_genStat(tr,         "MCParticles.generatorStatus");
  TTreeReaderArray<Int_t> mcpart_simStat(tr,         "MCParticles.simulatorStatus");
  TTreeReaderArray<Double_t> mcpart_m(tr,         "MCParticles.mass");
  TTreeReaderArray<Float_t> mcpart_psx(tr,       "MCParticles.momentum.x");
  TTreeReaderArray<Float_t> mcpart_psy(tr,       "MCParticles.momentum.y");
  TTreeReaderArray<Float_t> mcpart_psz(tr,       "MCParticles.momentum.z");


  // Reco tracks
  TTreeReaderArray<Int_t> recparts_type(tr,  "ReconstructedChargedParticles.type"); // needs to be made an int eventually in actual EE code
  TTreeReaderArray<Float_t> recparts_e(tr, "ReconstructedChargedParticles.energy");
  TTreeReaderArray<Float_t> recparts_p_x(tr, "ReconstructedChargedParticles.momentum.x");
  TTreeReaderArray<Float_t> recparts_p_y(tr, "ReconstructedChargedParticles.momentum.y");
  TTreeReaderArray<Float_t> recparts_p_z(tr, "ReconstructedChargedParticles.momentum.z");
  TTreeReaderArray<Int_t> recparts_PDG(tr,  "ReconstructedChargedParticles.PDG");
  TTreeReaderArray<Float_t> recparts_CHI2PID(tr,  "ReconstructedChargedParticles.goodnessOfPID");
  
  // RecoAssociations
  std::string assoc_branch_name = "ReconstructedChargedParticleAssociations";
  if(listOfBranches->FindObject(assoc_branch_name.c_str()) == nullptr)
    assoc_branch_name = "ReconstructedChargedParticlesAssociations"; // productions before 23.5
  TTreeReaderArray<UInt_t> assoc_simID(tr, (assoc_branch_name+".simID").c_str());
  TTreeReaderArray<UInt_t> assoc_recID(tr, (assoc_branch_name+".recID").c_str());
  TTreeReaderArray<Float_t> assoc_weight(tr, (assoc_branch_name+".weight").c_str());

  // calculate Q2 weights
  CalculateEventQ2Weights();

  // counters
  Long64_t numNoBeam, numEle, numNoEle, numNoHadrons, numProxMatched;
  numNoBeam = numEle = numNoEle = numNoHadrons = numProxMatched = 0;

  

  cout << "begin event loop..." << endl;
  tr.SetEntriesRange(1,maxEvents);
  do{
    if(tr.GetCurrentEntry()%10000==0) cout << tr.GetCurrentEntry() << " events..." << endl;

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();
    kin->vecHadrons.clear();
    kin->vecHadronsPIDs.clear();
    kinTrue->vecHadrons.clear();
    kinTrue->vecHadronsPIDs.clear();
    
    double maxP = 0;
    int genEleID = -1;
    bool foundBeamElectron = false;
    bool foundBeamIon = false;

    // Index maps for particle sets
    std::map <double,int> genidmap; // <pz, index_gen>
    std::map <int,int> mcidmap; // <index_mc, index_gen>
    std::map <int,int> recidmap; // <index, index_mc>  
    std::map <int,int> trackstatmap; // <index, genstatus>

    // Particles vectors
    // The index of the vectors correspond to their for loop idx
    std::vector<Particles> genpart;    // mcID --> igen
    std::vector<Particles> mcpart;     // mcID --> imc
    std::vector<Particles> recpart;  // mcID --> (imc of matching mcpart) or (-1 if no match is found)

    /*
      GenParticles loop
    */

    for(int igen=0; igen<hepmcp_PDG.GetSize(); igen++) {

      int pid_ = hepmcp_PDG[igen];

      double px_ = hepmcp_psx[igen];
      double py_ = hepmcp_psy[igen];
      double pz_ = hepmcp_psz[igen];
      double e_  = hepmcp_E[igen];
     
      double p_ = sqrt(pow(hepmcp_psx[igen],2) + pow(hepmcp_psy[igen],2) + pow(hepmcp_psz[igen],2));
      double mass_ = (fabs(pid_)==211)?constants::pimass:(fabs(pid_)==321)?constants::kmass:(fabs(pid_)==11)?constants::emass:(fabs(pid_)==13)?constants::mumass:(fabs(pid_)==2212)?constants::pmass:0.;

      // Add to genpart
      Particles part;

      part.pid=pid_;
      // part.charge = // TODO; not used yet
      part.mcID=igen;
      part.vecPart.SetPxPyPzE(px_,py_,pz_,e_);
      genpart.push_back(part);
      
      genidmap.insert({pz_,igen});

    }

    /*
      MCParticles loop
    */

    for(int imc=0; imc < mcpart_PDG.GetSize(); imc++){

      int pid_ = mcpart_PDG[imc];
      double px_ = mcpart_psx[imc];
      double py_ = mcpart_psy[imc];
      double pz_ = mcpart_psz[imc];
      double m_ = mcpart_m[imc];
      double e_ = sqrt(px_*px_+py_*py_+pz_*pz_+m_*m_);

      // Add to mcpart
      Particles part;

      part.pid=pid_;
      // part.charge = // TODO; not used yet
      part.mcID=imc;
      part.vecPart.SetPxPyPzE(px_,py_,pz_,e_);
      mcpart.push_back(part);
      
      int igen=-1;
      if(auto search = genidmap.find(pz_); search != genidmap.end())
	igen=search->second; //index of the GeneratedParticle
      
      mcidmap.insert({imc,igen});

    }

    /*
      ReconstructedParticles loop
      - Add all particles to the std::vector<> of particles
      - Look up associated MC particle
    */

   
      
    for(int irec=0; irec < recparts_PDG.GetSize(); irec++){

      int pid_ = recparts_PDG[irec];
      double px_ = recparts_p_x[irec];
      double py_ = recparts_p_y[irec];
      double pz_ = recparts_p_z[irec];
      double e_ = recparts_e[irec];
      
      // Add to recpart
      Particles part;

      part.pid=pid_;
      // part.charge = // TODO; not used yet
      part.vecPart.SetPxPyPzE(px_,py_,pz_,e_);

      double m_ = part.vecPart.M();

      /*
	Read through Associations to match particles
	By default, we assume no association, so mcID --> -1
	assoc_recID --> irec (index of the RecoParticle)
	assoc_simID --> imc (index of the MCParticle)
      */


      part.mcID=-1;
      for(int iassoc = 0 ; iassoc < assoc_simID.GetSize() ; iassoc++){
	int idx_recID = assoc_recID[iassoc]; 
	int idx_simID = assoc_simID[iassoc];
	if(irec==idx_recID){ // This track has an association
	  part.mcID=idx_simID;
	  break; // Only one association per particle
	}
      }

      recpart.push_back(part);
      recidmap.insert({irec,part.mcID}); 
      
    }



    /*
      With the GeneratedParticles, MCParticles, and ReconstructedParticles filled,
      we can begin to search for the beam particles and hadronic final state (HFS)
      This is done for both the Truth and Reconstructed Particles
    */

    /*
      Loop over MCParticles
    */

    for(auto mcpart_: mcpart){

      int imc = mcpart_.mcID;
      /* Beam particles have a MCParticles.generatorStatus of 4 */
      int genStat_ = mcpart_genStat[imc];
      if(mcpart_.pid==11 && genStat_ == 4){
	foundBeamElectron=true;
	kinTrue->vecEleBeam = mcpart_.vecPart;
      }
      else if(mcpart_.pid==2212 && genStat_ == 4){
	foundBeamIon=true;
	kinTrue->vecIonBeam = mcpart_.vecPart;
      }
      else if(genStat_==4){
	cout << "Warning...unknown beam particle with generatorStatus == 4 found...Continuing..." << endl;
      }

      /* Assume the scattered electron is the pid==11 final state particle with the most energy */
      if(mcpart_.pid==11 && genStat_ == 1 && mcpart_.vecPart.P() > maxP)
	{
	  maxP=mcpart_.vecPart.P();
	  kinTrue->vecElectron = mcpart_.vecPart;
	  genEleID = mcpart_.mcID; 
	}

      /*
	Only append MCParticles to the HFS if they are matched with a GeneratedParticle
       */
      
      else if(genStat_ == 1 && mcidmap[mcpart_.mcID]>-1){ 
	kinTrue->AddToHFS(mcpart_.vecPart);
      }
    
    }

    //check beam finding
    if(!foundBeamElectron || !foundBeamIon) { numNoBeam++; continue;};

    /*
      Loop over RecoParticles
    */

    int irec = 0;
    bool recEleFound=false;
    for(auto recpart_ : recpart){
      // If there is a matching MCParticle
      if(recidmap[irec]!=-1){
	// If the recidmap is linked to the genEleID (generated scattered electron), identify this reco particle as the electron
	if(recidmap[irec]==genEleID){	
	  recEleFound=true;
	  kin->vecElectron= recpart_.vecPart;
	}
	// Add the final state particle to the HFS
	kin->AddToHFS(recpart_.vecPart);

	// Add reconstructed particle and true info to HFSTree
	if( writeHFSTree ){
	  kin->AddToHFSTree(recpart_.vecPart, pid);	    
	}
      }
      irec++; // Increment to next particle
    }

    // Skip event if the reco scattered electron was missing
    if(recEleFound==false){
      numNoEle++;
      continue;
    }
    else
      numEle++;

    // subtrct electron from hadronic final state variables
    kin->SubtractElectronFromHFS();
    kinTrue->SubtractElectronFromHFS();

    // skip the event if there are no reconstructed particles (other than the
    // electron), otherwise hadronic recon methods will fail
    if(kin->countHadrons == 0) {
      numNoHadrons++;
      continue;
    };
   
    // calculate DIS kinematics
    if(!(kin->CalculateDIS(reconMethod))) continue; // reconstructed
    if(!(kinTrue->CalculateDIS(reconMethod))) continue; // generated (truth)

    // Get the weight for this event's Q2
    auto Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);
    
    // fill inclusive histograms, if only `inclusive` is included in output
    // (otherwise they will be filled in track and jet loops)
    if(includeOutputSet["inclusive_only"]) {
      auto wInclusive = Q2weightFactor * weightInclusive->GetWeight(*kinTrue);
      wInclusiveTotal += wInclusive;
      FillHistosInclusive(wInclusive);
    }

    /*
      Loop again over the reconstructed particles
      Calculate Hadron Kinematics
      Fill output data structures (Histos)
    */
    
    int ipart = 0;
    for(auto part : recpart){

      int pid_ = part.pid;
      int mcid_ = part.mcID;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      
      // First we start with the single particle tracks
      bool keep_reco_particle = false;
      auto kv = PIDtoFinalState.find(pid_);
      if(kv!=PIDtoFinalState.end()){
	finalStateID = kv->second;
	if(activeFinalStates.find(finalStateID)!=activeFinalStates.end()) keep_reco_particle=true;
      }
      // Second, we check the multiparticle tracks (ex: dihadrons)
      for (const auto& entry : PIDStoFinalState) {
	const std::vector<int>& key = entry.first;
	if (std::find(key.begin(), key.end(), pid_) != key.end()) {
	  finalStateID = entry.second;
	  if(activeFinalStates.find(finalStateID)!=activeFinalStates.end()) keep_reco_particle=true;
	}
      }
      
      // Now, we skip the particles which were not needed
      if(!keep_reco_particle) continue;

      // calculate reconstructed hadron kinematics
      kin->vecHadron = part.vecPart;
      kin->CalculateHadronKinematics();

      // if we are analyzing dihadrons, store the particle into a vector
      if(includeOutputSet["2h"]){
	kin->vecHadrons.push_back(part.vecPart);
	kin->vecHadronsPIDs.push_back(pid_);
      }
      // add selected single hadron FS to HFS tree
      if( writeHFSTree ){
	kin->AddTrackToHFSTree(part.vecPart, part.pid);
      }   
      
      
      // find the matching truth hadron using mcID, and calculate its kinematics
      if(mcid_ >= 0) {
	bool found_MC_match = false;
	for(auto imc : mcpart) {
	  if(mcid_ == imc.mcID) {
	    found_MC_match = true;
	    kinTrue->vecHadron = imc.vecPart;
	    // add tracks of interest for kinematic studies to HFSTree
	    if( writeHFSTree ){
	      kinTrue->AddTrackToHFSTree(imc.vecPart, imc.pid);				     
	    }
	    break;
	  }
	}
	if(includeOutputSet["2h"]){ // if analyzing dihadrons, append to list
	  if(found_MC_match)
	    kinTrue->vecHadrons.push_back(kinTrue->vecHadron);
	  else
	    kinTrue->vecHadrons.push_back(TLorentzVector());
	}
      }

      kinTrue->CalculateHadronKinematics();

      if(includeOutputSet["1h"]) {
        // fill single-hadron histograms in activated bins
        auto wTrack = Q2weightFactor * weightTrack->GetWeight(*kinTrue);
        wTrackTotal += wTrack;
        FillHistos1h(wTrack);
        FillHistosInclusive(wTrack);

	// fill simple tree
	// - not binned
	// - `IsActiveEvent()` is only true if at least one bin gets filled for this track
	if( writeSidisTree && HD->IsActiveEvent() ) ST->FillTree(wTrack);
      }
    } //hadron loop

    // Calculate dihadron kinematics
    if(includeOutputSet["2h"]) {
      // First hadron loop
      for(unsigned int i = 0 ; i < kin->vecHadrons.size(); i++){
	// Skip hadron if its PID isn't the first track wanted by user
	int pid_h1 = kin->vecHadronsPIDs.at(i);
	if(pid_h1==211) // *** Temporary code for only PiPlus PiMinus dihadrons ***
	  {
	    kin->vecDihadronA = kin->vecHadrons.at(i);
	    kinTrue->vecDihadronA = kinTrue->vecHadrons.at(i);
	  }
	else
	  continue;
	// Second hadron loop
	for(unsigned int j = 0 ; j < kin->vecHadrons.size(); j++){
	  // Skip hadron if its PID isn't the first track wanted by user
	  int pid_h2 = kin->vecHadronsPIDs.at(j);
	  if(pid_h2==-211) // *** Temporary code for only PiPlus PiMinus dihadrons ***
	    {
	      kin->vecDihadronB = kin->vecHadrons.at(j);
	      kinTrue->vecDihadronB = kinTrue->vecHadrons.at(j);
	    }
	  else
	    continue;

	  kin->CalculateDihadronKinematics();
	  kinTrue->CalculateDihadronKinematics();
	  
	  // Fill bins
	  auto wTrack = Q2weightFactor * weightTrack->GetWeight(*kinTrue);
	  FillHistos2h(wTrack); 
	  if( writeDiSidisTree && HD->IsActiveEvent() ) DiST->FillTree(wTrack);
	}	  
      }
    } // dihadron post loop
    
    if( writeHFSTree && kin->nHFS > 0) HFST->FillTree(Q2weightFactor);
    
    /*
      Loop again over the reconstructed particles
      Fill output data structures (ParticleTree)
    */

    // fill particle tree
    if ( writeParticleTree ) {
      int ipart = 0;
      for( auto part: recpart ){
	  Particles mcpart_;                  
	  int mcpart_idx=recidmap[ipart];     // Map idx to the matched MCParticle
	  int genStat_ = -1;                    // Default Generator Status of MCParticle is -1 (no match)
	  if(mcpart_idx>-1){                    // RecoParticle has an MCParticle match
	    mcpart_ = mcpart.at(mcpart_idx);        // Get MCParticle
	    int imc = mcpart_.mcID;          
	    genStat_ = mcpart_genStat[imc];         // Get Generator status of MCParticle
	    if(imc==genEleID)                       // If MCParticle was scattered electron, set status to 2
	      genStat_=2;
	  }
	  PT->FillTree(part.vecPart,      // Fill Tree
		       mcpart_.vecPart,
		       part.pid,
		       genStat_);
      } // particle loop
      ipart++;
    }
    
  } while(tr.Next());

  
  // finish execution
  Finish();

  // final printout
  cout << "Total number of scattered electrons found: " << numEle << endl;
  if(numNoEle>0)
    cerr << "WARNING: skipped " << numNoEle << " events which had no reconstructed scattered electron" << endl;
  if(numNoHadrons>0)
    cerr << "WARNING: skipped " << numNoHadrons << " events which had no reconstructed hadrons" << endl;
  if(numNoBeam>0)
    cerr << "WARNING: skipped " << numNoBeam << " events which had no beam particles" << endl;
  if(numProxMatched>0)
    cerr << "WARNING: " << numProxMatched << " recon. particles were proximity matched to truth (when mcID match failed)" << endl;
}
