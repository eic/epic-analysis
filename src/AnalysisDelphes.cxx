#include "AnalysisDelphes.h"

ClassImp(AnalysisDelphes)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
AnalysisDelphes::AnalysisDelphes(
  TString infileName_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_,
  TString outfilePrefix_
) : Analysis(
  infileName_,
  eleBeamEn_,
  ionBeamEn_,
  crossingAngle_,
  outfilePrefix_
) {
  // delphes-specific settings defaults
  /* ... none defined yet ... */
};

//=============================================
// perform the analysis
//=============================================
void AnalysisDelphes::Execute() {

  // setup
  Prepare();

  // read delphes tree
  TChain *chain = new TChain("Delphes");
  for(TString in : infiles) chain->Add(in);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  ENT = tr->GetEntries();

  // calculate cross section
  if(maxEvents>0) ENT = maxEvents; // limiter

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itPIDSystemsTrack(tr->UseBranch("PIDSystemsTrack"));

  // get counts in different Q2 cuts
  cout << "count events..." << endl;
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%10000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
    tr->ReadEntry(e);
    itParticle.Reset();
    maxElePtrue = 0;
    while(GenParticle *part = (GenParticle*) itParticle()){
      if(part->PID == 11 && part->Status == 1){
        elePtrue = part->PT * TMath::CosH(part->Eta);
        if(elePtrue > maxElePtrue){
          maxElePtrue = elePtrue;
          kinTrue->vecElectron.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              Kinematics::ElectronMass()
              );
        };
      };
    };
    kinTrue->CalculateDIS(reconMethod);
    Double_t Q2 = kinTrue->Q2;
    CountEvent(Q2, chain->GetTreeNumber());
  }

  cout << "count results:" << endl;
  for(Int_t idx=0; idx<inXsecs.size(); ++idx) {
    cout << "\tQ2 > " << inQ2mins[idx]
      << ": xs=" << inXsecs[idx] << ", n=" << inEntries[idx] << endl;
  }

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
    while(GenParticle *part = (GenParticle*) itParticle()){
      if(part->PID == 11 && part->Status == 1){
        elePtrue = part->PT * TMath::CosH(part->Eta);
        if(elePtrue > maxElePtrue){
          maxElePtrue = elePtrue;
          kinTrue->vecElectron.SetPtEtaPhiM(
              part->PT,
              part->Eta,
              part->Phi,
              Kinematics::ElectronMass()
              );
        };
      };
    };

    // get hadronic final state variables
    kin->GetHadronicFinalState(itTrack, itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);

    // calculate DIS kinematics
    kin->CalculateDIS(reconMethod); // reconstructed
    kinTrue->CalculateDIS(reconMethod); // generated (truth)
    Double_t Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, chain->GetTreeNumber());

    // get vector of jets
    // TODO: should this have an option for clustering method?
    kin->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
    
    // asymmetry injection
    //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
    //kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
    //kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

    // Get index of file that the event comes from.
    wTrack = Q2weightFactor * weight->GetWeight(*kinTrue);
    wTrackTotal += wTrack;


    // track loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    itTrack.Reset();
    while(Track *trk = (Track*) itTrack()) {
      //cout << e << " " << trk->PID << endl;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      pid = trk->PID;
      auto kv = PIDtoFinalState.find(pid);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      // get parent particle, to check if pion is from vector meson
      GenParticle *trkParticle = (GenParticle*)trk->Particle.GetObject();
      TObjArray *brParticle = (TObjArray*)itParticle.GetCollection();
      GenParticle *parentParticle = (GenParticle*)brParticle->At(trkParticle->M1);
      int parentPID = (parentParticle->PID); // TODO: this is not used yet...

      // calculate hadron kinematics
      kin->vecHadron.SetPtEtaPhiM(
          trk->PT,
          trk->Eta,
          trk->Phi,
          trk->Mass /* TODO: do we use track mass here ?? */
          );
      GenParticle* trkPart = (GenParticle*)trk->Particle.GetObject();
      kinTrue->vecHadron.SetPtEtaPhiM(
          trkPart->PT,
          trkPart->Eta,
          trkPart->Phi,
          trkPart->Mass /* TODO: do we use track mass here ?? */
          );

      kin->CalculateHadronKinematics();
      kinTrue->CalculateHadronKinematics();

      // fill track histograms in activated bins
      FillHistosTracks();

      // fill simple tree
      // - not binned
      // - `activeEvent` is only true if at least one bin gets filled for this track
      // - TODO [critical]: add a `finalState` cut (also needed in AnalysisDD4hep)
      if( writeSimpleTree && activeEvent ) ST->FillTree(wTrack);

    }; // end track loop


    // jet loop - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    finalStateID = "jet";
    if(activeFinalStates.find(finalStateID)!=activeFinalStates.end()) {

      #if INCCENTAURO == 1
      if(useBreitJets) kin->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
      #endif

      wJet = Q2weightFactor * weightJet->GetWeight(*kinTrue); // TODO: should we separate weights for breit and non-breit jets?
      wJetTotal += wJet;

      Int_t nJets;
      if(useBreitJets) nJets = kin->breitJetsRec.size();
      else      nJets = kin->jetsRec.size();

      for(int i = 0; i < kin->jetsRec.size(); i++){

        if(useBreitJets) {
          #if INCCENTAURO == 1
          jet = kin->breitJetsRec[i];
          kin->CalculateBreitJetKinematics(jet);
          #endif
        } else {
          jet = kin->jetsRec[i];
          kin->CalculateJetKinematics(jet);
        };

        // fill jet histograms in activated bins
        FillHistosJets();

      };
    }; // end jet loop

  };
  cout << "end event loop" << endl;
  // event loop end =========================================================


  // finish execution
  Finish();
};


// destructor
AnalysisDelphes::~AnalysisDelphes() {
};

