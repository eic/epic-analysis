// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks, Connor Pecar, Duane Byer, Matthew McEneaney, Brian Page

#include "AnalysisDelphes.h"

ClassImp(AnalysisDelphes)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
AnalysisDelphes::AnalysisDelphes(TString configFileName_, TString outfilePrefix_) :
  Analysis(configFileName_, outfilePrefix_)
{
  // delphes-specific settings defaults
  /* ... none defined yet ... */
};

// destructor
AnalysisDelphes::~AnalysisDelphes() { };


//=============================================
// perform the analysis
//=============================================
void AnalysisDelphes::Execute() {

  // setup
  Prepare();

  // read delphes tree
  auto chain = std::make_unique<TChain>("Delphes");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      // std::cout << "Adding " << infiles[idx][idxF] << " with " << inEntries[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str(), inEntries[idx][idxF]);
    }
  }
  chain->CanDeleteRefs();
  auto tr = std::make_unique<ExRootTreeReader>(chain.get());
  ENT = tr->GetEntries();
  if(maxEvents>0) ENT = std::min(maxEvents,ENT);

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itpfRICHTrack(tr->UseBranch("pfRICHTrack"));
  TObjArrayIter itDIRCepidTrack(tr->UseBranch("barrelDIRC_epidTrack"));
  TObjArrayIter itDIRChpidTrack(tr->UseBranch("barrelDIRC_hpidTrack"));
  TObjArrayIter itBTOFepidTrack(tr->UseBranch("BTOF_eTrack"));
  TObjArrayIter itBTOFhpidTrack(tr->UseBranch("BTOF_hTrack"));
  TObjArrayIter itdualRICHagTrack(tr->UseBranch("dualRICHagTrack"));
  TObjArrayIter itdualRICHcfTrack(tr->UseBranch("dualRICHcfTrack"));

  // calculate Q2 weights
  CalculateEventQ2Weights();

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%10000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
    tr->ReadEntry(e);

    // electron loop
    // - finds max-momentum electron
    itElectron.Reset();
    maxEleP = 0;
    while(Electron *ele = (Electron*) itElectron()) {
      eleP = ele->PT * TMath::CosH(ele->Eta);
      if(eleP>maxEleP) {
        maxEleP = eleP;
        kin->vecElectron.SetPtEtaPhiM(
            ele->PT,
            ele->Eta,
            ele->Phi,
            Kinematics::ElectronMass()
            );
      };
    };
    if(maxEleP<0.001) continue; // no scattered electron found

    // - repeat for truth electron
    itParticle.Reset();
    maxElePtrue = 0;
    bool found_elec = false;
    bool found_ion = false;
    while(GenParticle *part = (GenParticle*) itParticle()){
      if(part->PID == 11 && part->Status == 1){
        elePtrue = part->PT * TMath::CosH(part->Eta);
        if(elePtrue > maxElePtrue){
          maxElePtrue = elePtrue;
          kinTrue->vecElectron.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        };
      };
      if(part->PID == 11 && part->Status == 4){
        if(!found_elec){
          found_elec = true;
          kinTrue->vecEleBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        }else{
          ErrorPrint("ERROR: Found two beam electrons in one event");
        };
      };
      if(part->PID != 11 && part->Status == 4){
        if(!found_ion){
          found_ion = true;
          kinTrue->vecIonBeam.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              part->Mass
              );
        }else{
          ErrorPrint("ERROR: Found two beam ions in one event");
        };
      };
    };
    if(!found_elec){
      ErrorPrint("ERROR: Didn't find beam electron in event");
    };
    if(!found_ion){
      ErrorPrint("ERROR: Didn't find beam ion in event");
    }

    // get hadronic final state variables
    kin->GetHFS(
        itTrack,
        itEFlowTrack,
        itEFlowPhoton,
        itEFlowNeutralHadron,
        itpfRICHTrack,
        itDIRCepidTrack, itDIRChpidTrack,
        itBTOFepidTrack, itBTOFhpidTrack,
        itdualRICHagTrack, itdualRICHcfTrack
        );
    kinTrue->GetTrueHFS(itParticle);

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

    // track loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    itTrack.Reset();
    while(Track *trk = (Track*) itTrack()) {
      //cout << e << " " << trk->PID << endl;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      // pid = trk->PID; //NOTE: trk->PID is currently not smeared so it just returns the truth-level PID
      pid = kin->GetTrackPID( // get smeared PID
          trk,
          itpfRICHTrack,
          itDIRCepidTrack, itDIRChpidTrack,
          itBTOFepidTrack, itBTOFhpidTrack,
          itdualRICHagTrack, itdualRICHcfTrack
          );
      auto kv = PIDtoFinalState.find(pid);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      // get parent particle, to check if pion is from vector meson
      GenParticle *trkParticle = (GenParticle*)trk->Particle.GetObject();
      TObjArray *brParticle = (TObjArray*)itParticle.GetCollection();
      GenParticle *parentParticle = (GenParticle*)brParticle->At(trkParticle->M1);
      int parentPID = (parentParticle->PID); // TODO: this is not used yet...

      // calculate hadron kinematics
      kin->hadPID = pid;
      kin->vecHadron.SetPtEtaPhiM(
          trk->PT,
          trk->Eta,
          trk->Phi,
          trk->Mass /* TODO: do we use track mass here ?? */
          );
      GenParticle* trkPart = (GenParticle*)trk->Particle.GetObject();
      kinTrue->hadPID = pid;
      kinTrue->vecHadron.SetPtEtaPhiM(
          trkPart->PT,
          trkPart->Eta,
          trkPart->Phi,
          trkPart->Mass /* TODO: do we use track mass here ?? */
          );
      
      kin->CalculateHadronKinematics();
      kinTrue->CalculateHadronKinematics();

      // asymmetry injection
      //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      //kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      //kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin
  
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

      // tests
      //kin->ValidateHeadOnFrame();

    }; // end track loop


    // jet loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(includeOutputSet["jets"]) {

      // get vector of jets
      // TODO: should this have an option for clustering method?
      //kin->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
      kinJet->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle, jetAlg, jetRad, jetMin);

      finalStateID = "jet";

#ifdef INCLUDE_CENTAURO
      if(useBreitJets) kinJet->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
#endif

      auto wJet = Q2weightFactor * weightJet->GetWeight(*kinJetTrue); // TODO: should we separate weights for breit and non-breit jets?
      wJetTotal += wJet;

      Int_t nJets;
      if(useBreitJets) nJets = kinJet->breitJetsRec.size();
      else      nJets = kinJet->jetsRec.size();

      for(int i = 0; i < kinJet->jetsRec.size(); i++){

        if(useBreitJets) {
#ifdef INCLUDE_CENTAURO
          jet = kinJet->breitJetsRec[i];
          kinJet->CalculateBreitJetKinematics(jet);
#endif
        } else {
          jet = kinJet->jetsRec[i];
          kinJet->CalculateJetKinematics(jet);

	  // Match Reco Jet to Nearest Truth Jet (Specify DR matching limit between jets)
	  kinJet->CalculateJetResolution(jetMatchDR);
        };

        // fill jet histograms in activated bins
        FillHistosJets(wJet);
        FillHistosInclusive(wJet);

      };
    }; // end jet loop

  };
  cout << "end event loop" << endl;
  // event loop end =========================================================


  // finish execution
  Finish();
  //cout << "DEBUG PID in HFS: nSmeared=" << kin->countPIDsmeared << "  nNotSmeared=" << kin->countPIDtrue << endl;
};
