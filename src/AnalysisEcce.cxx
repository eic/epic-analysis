#include "AnalysisEcce.h"

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
AnalysisEcce::AnalysisEcce(TString configFileName_, TString outfilePrefix_) :
  Analysis(configFileName_, outfilePrefix_),
  trackSource(0) // default track source is "all tracks"
{ };

// destructor
AnalysisEcce::~AnalysisEcce() { };


//=============================================
// perform the analysis
//=============================================
void AnalysisEcce::Execute()
{
  // setup
  Prepare();

  // read EventEvaluator tree
  TChain *chain = new TChain("event_tree");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      // std::cout << "Adding " << infiles[idx][idxF] << " with " << inEntries[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str(), inEntries[idx][idxF]);
    }
  }

  TTreeReader tr(chain);

  // Truth

  TTreeReaderArray<Int_t> hepmcp_status(tr, "hepmcp_status");
  TTreeReaderArray<Int_t> hepmcp_PDG(tr,    "hepmcp_PDG");
  TTreeReaderArray<Float_t> hepmcp_E(tr,      "hepmcp_E");
  TTreeReaderArray<Float_t> hepmcp_psx(tr,    "hepmcp_px");
  TTreeReaderArray<Float_t> hepmcp_psy(tr,    "hepmcp_py");
  TTreeReaderArray<Float_t> hepmcp_psz(tr,    "hepmcp_pz");

  TTreeReaderArray<Int_t> hepmcp_BCID(tr, "hepmcp_BCID");
  TTreeReaderArray<Int_t> hepmcp_m1(tr, "hepmcp_m1"); 
  TTreeReaderArray<Int_t> hepmcp_m2(tr, "hepmcp_m2"); 
  

  // All true particles (including secondaries, etc)
  TTreeReaderArray<Int_t> mcpart_ID(tr,        "mcpart_ID");
  TTreeReaderArray<Int_t> mcpart_ID_parent(tr, "mcpart_ID_parent"); 
  TTreeReaderArray<Int_t> mcpart_PDG(tr,       "mcpart_PDG");
  TTreeReaderArray<Float_t> mcpart_E(tr,         "mcpart_E");
  TTreeReaderArray<Float_t> mcpart_psx(tr,       "mcpart_px");
  TTreeReaderArray<Float_t> mcpart_psy(tr,       "mcpart_py");
  TTreeReaderArray<Float_t> mcpart_psz(tr,       "mcpart_pz");
  TTreeReaderArray<Int_t> mcpart_BCID(tr,      "mcpart_BCID");


  // Reco tracks
  TTreeReaderArray<Float_t> tracks_id(tr,  "tracks_ID"); // needs to be made an int eventually in actual EE code
  TTreeReaderArray<Float_t> tracks_p_x(tr, "tracks_px");
  TTreeReaderArray<Float_t> tracks_p_y(tr, "tracks_py");
  TTreeReaderArray<Float_t> tracks_p_z(tr,  "tracks_pz");
  TTreeReaderArray<Float_t> tracks_trueID(tr,  "tracks_trueID");
  TTreeReaderArray<UShort_t> tracks_source(tr,  "tracks_source");
   // TTreeReaderArray<Short_t> tracks_charge(tr,  "tracks_charge");


  // calculate Q2 weights
  CalculateEventQ2Weights();

  // counters
  Long64_t numNoBeam, numEle, numNoEle, numNoHadrons, numProxMatched;
  numNoBeam = numEle = numNoEle = numNoHadrons = numProxMatched = 0;

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  tr.SetEntriesRange(1,maxEvents);
  do {
    if(tr.GetCurrentEntry()%10000==0) cout << tr.GetCurrentEntry() << " events..." << endl;
  
    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();




    // a few maps needed to get the associated info between tracks, true particles, etec

    std::map <int,int> mcidmap;
    std::map <int,int> mcbcidmap;
    std::map <int,int> mcbcidmap2;    
    std::map <int,int> trackmap; // mapping of all the tracks
        
    for (int imc =0; imc < mcpart_ID.GetSize(); imc++){
      if (mcpart_E[imc]<0.1) continue;       
      int id = (int)mcpart_ID[imc];
      int bcid = (int)mcpart_BCID[imc];
      mcidmap.insert({id,imc});
      mcbcidmap.insert({bcid,imc});       
    }    
    
    for(int itrack=0; itrack<tracks_id.GetSize(); itrack++) {
      if(tracks_source[itrack]!=trackSource) continue;
      trackmap.insert({tracks_trueID[itrack],itrack});
    }
    




    // generated truth loop
    /* - add truth particle to `mcpart`
     * - add to hadronic final state sums (momentum, sigma, etc.)
     * - find scattered electron
     * - find beam particles
     */
    std::vector<ParticlesEE> mcpart;
    double maxPt = 0;
    int genEleID = -1;
    bool foundBeamElectron = false;
    bool foundBeamIon = false;
    int genEleBCID = -1;

    for(int imc=0; imc<hepmcp_PDG.GetSize(); imc++) {
      
      int pid_ = hepmcp_PDG[imc];
            
      int genStatus_ = hepmcp_status[imc]; // genStatus 4: beam particle,  1: final state
      
      double px_ = hepmcp_psx[imc];
      double py_ = hepmcp_psy[imc];
      double pz_ = hepmcp_psz[imc];
      double e_  = hepmcp_E[imc];
      
      double p_ = sqrt(pow(hepmcp_psx[imc],2) + pow(hepmcp_psy[imc],2) + pow(hepmcp_psz[imc],2));
      double pt_ = sqrt(pow(hepmcp_psx[imc],2) + pow(hepmcp_psy[imc],2));
      double mass_ = (fabs(pid_)==211)?pimass:(fabs(pid_)==321)?kmass:(fabs(pid_)==11)?emass:(fabs(pid_)==13)?mumass:(fabs(pid_)==2212)?pmass:0.;      

      // add to `mcpart`
      ParticlesEE part;

      //cout << genStatus_ << " " << pid_ << " " << part.mcID << " " << hepmcp_BCID[imc] << " " << hepmcp_m1[imc] << " " << hepmcp_m2[imc] << " "
      //      << px_ << " " << py_ << " " << pz_ << " " << e_ << endl;

      
      if(genStatus_ == 1) { // final state
	
	int imcpart = -1;// matched truthtrack
	auto search = mcbcidmap.find((int)(hepmcp_BCID[imc]));
	if (search != mcbcidmap.end()) {
	  imcpart = search->second;
	}
	
	
	if (imcpart >-1){
	  px_ = mcpart_psx[imcpart];
	  py_ = mcpart_psy[imcpart];
	  pz_ = mcpart_psz[imcpart];
	  e_  = mcpart_E[imcpart];	  
	  p_ = sqrt(pow(mcpart_psx[imcpart],2) + pow(mcpart_psy[imcpart],2) + pow(mcpart_psz[imcpart],2));
	  pt_ = sqrt(pow(mcpart_psx[imcpart],2) + pow(mcpart_psy[imcpart],2));
	  part.mcID = mcpart_ID[imcpart];
	}
	  else
	    part.mcID = -1;
	//cout << genStatus_ << " " << pid_ << " " << part.mcID << " " << hepmcp_BCID[imc] << " " << hepmcp_m1[imc] << " " << hepmcp_m2[imc] << " "
	//    << px_ << " " << py_ << " " << pz_ << " " << e_ << endl;
	
	  part.pid = pid_;
	  part.vecPart.SetPxPyPzE(px_, py_, pz_, e_);
        
	  mcpart.push_back(part);
	
	  // add to hadronic final state sums
	  kinTrue->AddToHFS(part.vecPart);


	  // identify scattered electron by max momentum
	  if(pid_ == 11) {
	    if(pt_ > maxPt) {
	      maxPt = pt_;
	      kinTrue->vecElectron.SetPxPyPzE(px_, py_, pz_, e_);
	      genEleID = part.mcID; //mcpart_ID[imc];
	      genEleBCID = hepmcp_BCID[imc];
	      //	      cout  << "\t\t\t found scattered electron  " << Form(" %6.2f %6.2f %6.2f %6.2f %6.2f  %5.3f %6.2f %6.2f id %3d\n",px_,py_,pz_, sqrt(p_*p_ + mass_*mass_),p_,mass_,hepmcp_E[imc],mcpart_E[imcpart],genEleID);
	    }
	  }
      }
      
      else if(genStatus_ == 4) { // beam particles
        if(pid_ == 11) { // electron beam
          if(!foundBeamElectron) {
            foundBeamElectron = true;
            kinTrue->vecEleBeam.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
	    //	    cout  << "\t\t\t found beam electron  " << Form(" %4.2f %4.2f %4.2f \n",px_,py_,pz_);
          }
          else { ErrorPrint("ERROR: Found two beam electrons in one event"); }
        }
        else { // ion beam
          if(!foundBeamIon) {
            foundBeamIon = true;
            kinTrue->vecIonBeam.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
	    //	    cout  << "\t\t\t found beam ion  " << Form(" %4.2f %4.2f %4.2f \n",px_,py_,pz_);
          }
          else { ErrorPrint("ERROR: Found two beam ions in one event"); }
        }
      }
    } // end truth loop
    
    for(int imc=0; imc<hepmcp_PDG.GetSize(); imc++) {
      int pid_ = hepmcp_PDG[imc];     
      int genStatus_ = hepmcp_status[imc]; // genStatus 4: beam particle,  1: final state                                                                            						

      double px_ = hepmcp_psx[imc];
      double py_ = hepmcp_psy[imc];
      double pz_ = hepmcp_psz[imc];
      double e_  = hepmcp_E[imc];
      
      double p_ = sqrt(pow(hepmcp_psx[imc],2) + pow(hepmcp_psy[imc],2) + pow(hepmcp_psz[imc],2));
      double mass_ = (fabs(pid_)==211)?pimass:(fabs(pid_)==321)?kmass:(fabs(pid_)==11)?emass:(fabs(pid_)==13)?mumass:(fabs(pid_)==2212)?pmass:0.;      
      // add to `mcpart`                                                                                                                                                                                                                      
	   int mother = hepmcp_m1[imc];
      if(pid_ == 22 || pid_ == 23){
	if(genStatus_ == 1){
	  if( mother == genEleBCID ){
	    TLorentzVector FSRmom(px_, py_, pz_, e_);
	    TLorentzVector eleCorr = kinTrue->vecElectron + FSRmom;
	    kinTrue->vecElectron.SetPxPyPzE(eleCorr.Px(), eleCorr.Py(), eleCorr.Pz(), eleCorr.E());
	  }
	}
      }
      
      
    }
    // check beam finding
    if(!foundBeamElectron || !foundBeamIon) { numNoBeam++; continue; };


    // reconstructed particles loop
    /* - add reconstructed particle to `recopart`
     * - find the scattered electron
     *
     */
    std::vector<ParticlesEE> recopart;
    int recEleFound = 0;
    for(int ireco=0; ireco<tracks_id.GetSize(); ireco++) {

      // skip track if not from the user-specified source `trackSource`
      if(tracks_source[ireco]!=trackSource) continue;

      int pid_ = 0; //tracks_pid[ireco];
      double reco_mass = 0.;
      int imc = -1;// matched truthtrack
      auto search = mcidmap.find((int)(tracks_trueID[ireco]));
      if (search != mcidmap.end()) {
	imc = search->second;
	pid_ = (int)(mcpart_PDG[imc]);
      }
      // later also use the likelihoods instead for pid
      //      cout  << "\t\t track " << ireco << " PDG  " << pid_ << " mcpart imc " << imc << endl;
      if(pid_ == 0) continue; // pid==0: reconstructed tracks with no matching truth pid

      // add reconstructed particle `part` to `recopart`
      ParticlesEE part;
      part.pid = pid_;
      part.mcID = tracks_trueID[ireco];
      //      part.charge = tracks_charge[ireco];
      part.charge = (pid_ == 211 || pid_ == 321 || pid_ == 2212 || pid_ == -11 || pid_ == -13)?1:(pid_ == -211 || pid_ == -321 || pid_ == -2212 || pid_ == 11 || pid_ == 13)?-1:0;

      double reco_px = tracks_p_x[ireco];
      double reco_py = tracks_p_y[ireco];
      double reco_pz = tracks_p_z[ireco];
      reco_mass = (fabs(pid_)==211)?pimass:(fabs(pid_)==321)?kmass:(fabs(pid_)==11)?emass:(fabs(pid_)==13)?mumass:(fabs(pid_)==2212)?pmass:0.;
      
      double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
      double reco_E = sqrt(reco_p*reco_p + reco_mass * reco_mass);

      part.vecPart.SetPxPyPzE(reco_px, reco_py, reco_pz, sqrt(reco_p*reco_p + reco_mass*reco_mass));


      //cout  << "\t\t\t track  " << Form(" %4.2f %4.2f %4.2f true id %4d imc %3d mcid %3d \n",reco_px,reco_py,reco_pz,tracks_trueID[ireco],imc,part.mcID);
      
      // add to `recopart` and hadronic final state sums only if there is a matching truth particle
      if(part.mcID > 0 && part.mcID != genEleID) { 
	if(imc>-1) {
	  //  cout  << "\t\t\t add  to hadfs  \n" ;
	  recopart.push_back(part);
	  kin->AddToHFS(part.vecPart);
	  kin->hfspid[kin->nHFS - 1] = pid_;
	}
      }

      // find scattered electron, by matching to truth // TODO: not realistic... is there an upstream electron finder?
      if(pid_ == 11 && part.mcID == genEleID) {
        recEleFound++;
	//	cout  << "\t\t\t reco electron " << Form(" %6.2f %6.2f %6.2f %6.2f  id %3d\n",reco_px,reco_py,reco_pz,sqrt(reco_p*reco_p + reco_mass*reco_mass),genEleID);
        kin->vecElectron.SetPxPyPzE(reco_px, reco_py, reco_pz, sqrt(reco_p*reco_p + reco_mass*reco_mass));
      }

    } // end reco loop

    // skip the event if the scattered electron is not found and we need it
    if(recEleFound < 1) {
      numNoEle++;
      continue; // TODO: only need to skip if we are using a recon method that needs it (`if reconMethod==...`)
    }
    else if(recEleFound>1) ErrorPrint("WARNING: found more than 1 reconstructed scattered electron in an event");
    else numEle++;

    // subtract electron from hadronic final state variables
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

    //fill HFS tree here?
    if( writeHFSTree && kin->nHFS > 0) HFST->FillTree(1.0);
    
    
    // loop over reconstructed particles again
    /* - calculate hadron kinematics
     * - fill output data structures (Histos, SimpleTree, etc.)
     */
    for(auto part : recopart) {
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
      if(mcid_ > 0) {
        for(auto imc : mcpart) {
          if(mcid_ == imc.mcID) {
            kinTrue->vecHadron = imc.vecPart;
            break;
          }
        }
      }
      /* // deprecated, since existence of truth match is checked earlier; in practice prox matching was never called
	 else {
	 // give it another shot: proximity matching
	 double mineta = 4.0;
	 numProxMatched++;
	 for(int imc=0; imc<(int)mcpart.size(); imc++) {
	 if(pid_ == mcpart[imc].pid) {
	 double deta = abs(kin->vecHadron.Eta() - mcpart[imc].vecPart.Eta());
            if(deta < mineta) {
              mineta = deta;
              kinTrue->vecHadron = mcpart[imc].vecPart;
            }
          }
        }
      }
      */
      kinTrue->CalculateHadronKinematics();

      // asymmetry injection
      // kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      // kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      // kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

      // weighting
      Double_t Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);
      wTrack = Q2weightFactor * weight->GetWeight(*kinTrue);
      wTrackTotal += wTrack;

      if(includeOutputSet["1h"]) {
        // fill track histograms in activated bins
        FillHistosTracks();

        // fill simple tree
        // - not binned
        // - `IsActiveEvent()` is only true if at least one bin gets filled for this track
        if( writeSimpleTree && HD->IsActiveEvent() ) ST->FillTree(wTrack);
      }

    }//hadron loop

  } while(tr.Next()); // tree reader loop
  cout << "end event loop" << endl;
  // event loop end =========================================================

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
