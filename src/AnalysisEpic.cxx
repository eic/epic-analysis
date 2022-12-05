// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Gregory Matousek, Christopher Dilks

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
  TChain *chain = new TChain("events");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      // std::cout << "Adding " << infiles[idx][idxF] << " with " << inEntries[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str(), inEntries[idx][idxF]);
    }
  }

  TTreeReader tr(chain);

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
  TTreeReaderArray<UInt_t> assoc_simID(tr, "ReconstructedChargedParticlesAssociations.simID");
  TTreeReaderArray<UInt_t> assoc_recID(tr, "ReconstructedChargedParticlesAssociations.recID");
  TTreeReaderArray<Float_t> assoc_weight(tr, "ReconstructedChargedParticlesAssociations.weight");

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
      part.charge = (pid_ == 211 || pid_ == 321 || pid_ == 2212 || pid_ == -11 || pid_ == -13)?1:(pid_ == -211 || pid_ == -321 || pid_ == -2212 || pid_ == 11 || pid_ == 13)?-1:0;
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
      part.charge = (pid_ == 211 || pid_ == 321 || pid_ == 2212 || pid_ == -11 || pid_ == -13)?1:(pid_ == -211 || pid_ == -321 || pid_ == -2212 || pid_ == 11 || pid_ == 13)?-1:0;
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
      - Identify the 
      - Identify closest matching MCParticle in theta,phi,E space

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
      part.charge = (pid_ == 211 || pid_ == 321 || pid_ == 2212 || pid_ == -11 || pid_ == -13)?1:(pid_ == -211 || pid_ == -321 || pid_ == -2212 || pid_ == 11 || pid_ == 13)?-1:0;
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
      auto kv = PIDtoFinalState.find(pid_);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      // calculate reconstructed hadron kinematics
      kin->vecHadron = part.vecPart;
      kin->CalculateHadronKinematics();

      // find the matching truth hadron using mcID, and calculate its kinematics
      if(mcid_ >= 0) {
	for(auto imc : mcpart) {
	  if(mcid_ == imc.mcID) {
	    kinTrue->vecHadron = imc.vecPart;
	    break;
	  }
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
	if( writeSimpleTree && HD->IsActiveEvent() ) ST->FillTree(wTrack);
      }
    } //hadron loop
    
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
