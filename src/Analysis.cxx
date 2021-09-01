#include "Analysis.h"

ClassImp(Analysis)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
Analysis::Analysis(
  TString infileName_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_,
  TString outfilePrefix_
)
  : infileName(infileName_)
  , eleBeamEn(eleBeamEn_)
  , ionBeamEn(ionBeamEn_)
  , crossingAngle(crossingAngle_)
  , outfilePrefix(outfilePrefix_)
{
  // build list of variables available for binning (paired with titles)
  // - availableBinSchemes is a map from variable name to variable title
  // - try to avoid using underscores in the variable name (they are okay in the title)
  availableBinSchemes.insert(std::pair<TString,TString>("pt","p_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("z","z"));
  availableBinSchemes.insert(std::pair<TString,TString>("x","x"));
  availableBinSchemes.insert(std::pair<TString,TString>("q","Q"));
  availableBinSchemes.insert(std::pair<TString,TString>("y","y"));
  availableBinSchemes.insert(std::pair<TString,TString>("ptJet", "jet p_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("zJet", "jet z"));
  availableBinSchemes.insert(std::pair<TString,TString>("finalState","finalState"));
  availableBinSchemes.insert(std::pair<TString,TString>("recMethod","recMethod"));

  /* TODO: re-enable
  // define final state bins (e.g., tracks or jets)
  AddBinScheme("finalState");
  AddFinalState("pipTrack","#pi^{+} tracks", 211);
  //AddFinalState("pimTrack","#pi^{-} tracks",-211);

  // define reconstruction method bins
  AddBinScheme("recMethod");
  AddRecMethod("Ele", "electron method");
  AddRecMethod("DA", "DA method");
  AddRecMethod("JB", "JB method");
  */

  /*
  // initialize diagonalizer settings
  // TODO: generalized diagonalizer
  diagonalPtXZ = false;
  diagonalXZQ = false;
  */

  // initializations
  writeSimpleTree = false;
  maxEvents = 0;
};


//=============================================
// perform the analysis
//=============================================
void Analysis::Execute() {

  cout << "-- running analysis of " << infileName << endl;

  // define output file
  outfileName = infileName;
  outfileName(TRegexp("^.*/")) = ""; // remove path
  outfileName(TRegexp("\\*")) = ""; // remove asterisk wildcard
  if(outfilePrefix!="") outfilePrefix+=".";
  outfileName = "out/"+outfilePrefix+outfileName;
  cout << "-- output file: " << outfileName << endl;
  outFile = new TFile(outfileName,"RECREATE");

  // instantiate objects
  kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);
  kinTrue = new Kinematics(eleBeamEn, ionBeamEn, crossingAngle);
  ST = new SimpleTree("tree",kin);
  weight = new WeightsUniform();
  weightJet = new WeightsUniform();
  HD = new HistosDAG();


  // read delphes tree
  TChain *chain = new TChain("Delphes");
  chain->Add(infileName);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  Long64_t ENT = tr->GetEntries();

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itPIDSystemsTrack(tr->UseBranch("PIDSystemsTrack"));

  // number of bins
  /*
  const Int_t NptBins = BinScheme("pt")->GetNumBins();
  const Int_t NxBins = BinScheme("x")->GetNumBins();
  const Int_t NzBins = BinScheme("z")->GetNumBins();
  const Int_t NqBins = BinScheme("q")->GetNumBins();
  const Int_t NyBins = BinScheme("y")->GetNumBins();
  const Int_t NfinalStateBins = BinScheme("finalState")->GetNumBins();
  const Int_t NptjetBins = BinScheme("ptJet")->GetNumBins();
  const Int_t NzjetBins = BinScheme("zJet")->GetNumBins();
  const Int_t NrecMethodBins = BinScheme("recMethod")->GetNumBins();
  */

  // sets of histogram sets
  // - `histSet` is a data structure for storing and organizing pointers to
  //   sets of histograms (`Histos` objects)
  // - `histSet*List` are used as temporary lists of relevant `Histos` pointers
  // - TODO: if we add one more dimension, 7D array will probably break; need
  //         better data structure

  HistosDAG *HD = new HistosDAG();
  HD->Build(binSchemes);

  /*
  Histos *histSet[NptBins][NxBins][NzBins][NqBins][NyBins][NfinalStateBins];
  Histos *histSetJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins];
  Histos *histSetBreitJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins][NrecMethodBins];
  */


  // HISTOGRAMS ================================================
  HD->ForEach([this](Histos *HS){
    // -- DIS kinematics
    HS->DefineHist2D("Q2vsX","x","Q^{2}","","GeV^{2}",
        NBINS,1e-3,1,
        NBINS,1,100,
        true,true
        );
    HS->DefineHist1D("Q","Q","GeV",NBINS,1.0,11.0,true,true);
    HS->DefineHist1D("x","x","",NBINS,1e-3,1.0,true,true);
    HS->DefineHist1D("y","y","",NBINS,1e-5,1,true);
    HS->DefineHist1D("W","W","GeV",NBINS,0,15);
    // -- DIS kinematics resolution
    HS->DefineHist1D("xRes","x - x_{true}","", NBINS, -1, 1);
    // -- hadron 4-momentum
    HS->DefineHist1D("pLab","p_{lab}","GeV",NBINS,0,10);
    HS->DefineHist1D("pTlab","p_{T}^{lab}","GeV",NBINS,1e-2,3,true);
    HS->DefineHist1D("etaLab","#eta_{lab}","",NBINS,-5,5);
    HS->DefineHist1D("phiLab","#phi_{lab}","",NBINS,-TMath::Pi(),TMath::Pi());
    // -- hadron kinematics
    HS->DefineHist1D("z","z","",NBINS,0,1);
    HS->DefineHist1D("pT","p_{T}","GeV",NBINS,1e-2,3,true);
    HS->DefineHist1D("qT","q_{T}","GeV",NBINS,1e-2,5,true);
    HS->DefineHist1D("qTq","q_{T}/Q","",NBINS,1e-2,3,true);
    HS->DefineHist1D("mX","m_{X}","GeV",NBINS,0,20);
    HS->DefineHist1D("phiH","#phi_{h}","",NBINS,-TMath::Pi(),TMath::Pi());
    HS->DefineHist1D("phiS","#phi_{S}","",NBINS,-TMath::Pi(),TMath::Pi());
    // -- cross sections
    //HS->DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
    HS->DefineHist1D("Q_xsec","Q","GeV",10,1.0,10.0,true,true); // log
    HS->Hist("Q_xsec")->SetMinimum(1e-10);
    // ===========================================================
  });
  HD->ExecuteAndClearOps();


  /*
  cout << "Define jet histograms..." << endl;
  for(int bpt=0; bpt<NptjetBins; bpt++) { // - loop over jet pT bins
    for(int bz=0; bz<NzjetBins; bz++){
      for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q bins
          for(int by=0; by<NyBins; by++) { // - loop over y bins

            histosN = this->GetHistosNameJets(bpt, bz, bx, bq, by);
            histosT = this->GetHistosTitleJets(bpt, bz, bx, bq, by);

            histSetJets[bpt][bz][bx][bq][by] = new Histos(histosN,histosT);
            HS = histSetJets[bpt][bz][bx][bq][by]; // shorthand pointer

            // jet kinematics plots
            HS->DefineHist1D("pT_jet","p_{T}","GeV", NBINS, 1e-2, 50);
            HS->DefineHist1D("mT_jet","m_{T}","GeV", NBINS, 1e-2, 20);
            HS->DefineHist1D("zJet","z","GeV", NBINS,0, 1);
            HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
            HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0);
            HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
            HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);
            // store cut definitions with histogram sets, then add histogram sets full list
            HS->AddCutDef(BinScheme("ptJet")->Cut(bpt));
            HS->AddCutDef(BinScheme("zJet")->Cut(bz));
            HS->AddCutDef(BinScheme("x")->Cut(bx));
            HS->AddCutDef(BinScheme("q")->Cut(bq));
            HS->AddCutDef(BinScheme("y")->Cut(by));
            histSetListJets.push_back(histSetJets[bpt][bz][bx][bq][by]);
          };	
        };
      };
    };
  };
  #if INCCENTAURO == 1
  for(int bpt=0; bpt<NptjetBins; bpt++) { // - loop over jet pT bins
    for(int bz=0; bz<NzjetBins; bz++){
      for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q bins
          for(int by=0; by<NyBins; by++) { // - loop over y bins
            for(int brec=0; brec<NrecMethodBins; brec++){
              histosN = this->GetHistosNameBreitJets(bpt, bz, bx, bq, by, brec);
              histosT = this->GetHistosTitleBreitJets(bpt, bz, bx, bq, by, brec);
              histSetBreitJets[bpt][bz][bx][bq][by][brec] = new Histos(histosN,histosT);
              HS = histSetBreitJets[bpt][bz][bx][bq][by][brec]; // shorthand pointer

              // jet kinematics plots
              HS->DefineHist1D("pT_jet","p_{T}","GeV", NBINS, 1e-2, 20);
              HS->DefineHist1D("mT_jet","m_{T}","GeV", NBINS, 1e-2, 20);
              HS->DefineHist1D("zJet","z","GeV", NBINS,0, 1);
              HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
              HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0);
              HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
              HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);

              // store cut definitions with histogram sets, then add histogram sets full list
              HS->AddCutDef(BinScheme("ptJet")->Cut(bpt));
              HS->AddCutDef(BinScheme("zJet")->Cut(bz));
              HS->AddCutDef(BinScheme("x")->Cut(bx));
              HS->AddCutDef(BinScheme("q")->Cut(bq));
              HS->AddCutDef(BinScheme("y")->Cut(by));
              HS->AddCutDef(BinScheme("recMethod")->Cut(brec));
              histSetListBreitJets.push_back(histSetBreitJets[bpt][bz][bx][bq][by][brec]);
            };
          };
        };
      };
    };
  };
  #endif
  */

  // get cross section and number of events
  // - cross sections are hard-coded, coped from pythia output
  Int_t eleBeamEnInt = (Int_t) eleBeamEn;
  Int_t ionBeamEnInt = (Int_t) ionBeamEn;
  Double_t xsecTot; // [nb]
  if     (eleBeamEnInt==5  && ionBeamEnInt==41 ) xsecTot=297.9259;
  else if(eleBeamEnInt==18 && ionBeamEnInt==275) xsecTot=700.0; // TODO; this is approximate
  else {
    cerr << "WARNING: unknown cross section; integrated lumi will be wrong" << endl;
    xsecTot=1;
  };
  Long64_t numGen = tr->GetEntries();
  TString sep = "--------------------------------------------";
  cout << sep << endl;
  cout << "assumed total cross section: " << xsecTot << " nb" << endl;
  cout << "number of generated events:  " << numGen << endl;

  // initialize total weights
  Double_t wTotal = 0.;
  Double_t wJetTotal = 0.;

  // vars
  Double_t eleP,maxEleP;
  int pid,bFinalState;
  Double_t elePtrue, maxElePtrue;

  // event loop =========================================================
  if(maxEvents>0) ENT = maxEvents; // limiter
  cout << "begin event loop..." << endl;
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%100000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
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

    maxElePtrue = 0;
    while(GenParticle *part = (GenParticle*) itParticle()){
      if(part->PID == 11 && part->Status == 1){
        elePtrue = part->PT;
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
    kin->CalculateDISbyElectron(); // reconstructed
    kinTrue->CalculateDISbyElectron(); // generated (truth)

    // get vector of jets
    // TODO: should this have an option for clustering method?
    kin->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);

    // track loop
    itTrack.Reset();
    while(Track *trk = (Track*) itTrack()) {
      //cout << e << " " << trk->PID << endl;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next
      pid = trk->PID;
      auto kv = PIDtoEnum.find(pid);
      if(kv!=PIDtoEnum.end()) bFinalState = kv->second;
      else continue;

      // get parent particle, to check if pion is from vector meson
      GenParticle *trkParticle = (GenParticle*)trk->Particle.GetObject();
      TObjArray *brParticle = (TObjArray*)itParticle.GetCollection();
      GenParticle *parentParticle = (GenParticle*)brParticle->At(trkParticle->M1);
      int parentPID = (parentParticle->PID);


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

      Double_t w = weight->GetWeight(*kin);
      wTotal += w;

      // APPLY MAIN CUTS
      if(kin->CutFull()) {

        // map varNames to values
        valueMap.clear();
        valueMap.insert(std::pair<TString,Double_t>(  "pt",  kin->pT  ));
        valueMap.insert(std::pair<TString,Double_t>(  "x",   kin->x   ));
        valueMap.insert(std::pair<TString,Double_t>(  "z",   kin->z   ));
        valueMap.insert(std::pair<TString,Double_t>(  "q2",  kin->Q2  ));
        valueMap.insert(std::pair<TString,Double_t>(  "y",   kin->y   ));

        // check which bins the event falls in
        Bool_t activeEvent = false;
        HD->TraverseBreadth([this,&activeEvent](Node *N){
          if(N->GetNodeType()==NT::bin) {
            auto val = valueMap.at(N->GetVarName());
            Bool_t active = N->GetCut()->CheckCut(val);
            if(active) activeEvent=true;
            N->SetActiveState(active);
          };
        });

        // fill histograms
        HD->ForEach([this,&w](Histos *H){
          // DIS kinematics
          dynamic_cast<TH2*>(H->Hist("Q2vsX"))->Fill(kin->x,kin->Q2,w);
          H->Hist("Q")->Fill(TMath::Sqrt(kin->Q2),w);
          H->Hist("x")->Fill(kin->x,w);
          H->Hist("W")->Fill(kin->W,w);
          H->Hist("y")->Fill(kin->y,w);
          // hadron 4-momentum
          H->Hist("pLab")->Fill(kin->pLab,w);
          H->Hist("pTlab")->Fill(kin->pTlab,w);
          H->Hist("etaLab")->Fill(kin->etaLab,w);
          H->Hist("phiLab")->Fill(kin->phiLab,w);
          // hadron kinematics
          H->Hist("z")->Fill(kin->z,w);
          H->Hist("pT")->Fill(kin->pT,w);
          H->Hist("qT")->Fill(kin->qT,w);
          H->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2),w);
          H->Hist("mX")->Fill(kin->mX,w);
          H->Hist("phiH")->Fill(kin->phiH,w);
          H->Hist("phiS")->Fill(kin->phiS,w);
          // cross sections (divide by lumi after all events processed)
          H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),w);
          // DIS kinematics resolution
          H->Hist("xRes")->Fill(kin->x - kinTrue->x,w);
        });
        HD->ExecuteOps(true); // save time and don't ClearOps (next loop will overwrite lambda)

        // fill simple tree (not binned)
        // TODO: consider adding a `finalState` cut
        if( writeSimpleTree && activeEvent ) ST->FillTree(w);

      };
    };

    // jet loop
    /*
    if(kin->CutDIS()){

      Double_t wJet = weightJet->GetWeight(*kin);
      wJetTotal += wJet;

      for(int i = 0; i < kin->jetsRec.size(); i++){
        PseudoJet jet = kin->jetsRec[i];
        kin->CalculateJetKinematics(jet);

        // following same procedure as in track loop	
        CheckBins( BinScheme("ptJet"), v_pt, kin->pTjet );
        CheckBins( BinScheme("zJet"), v_z, kin->zjet );
        CheckBins( BinScheme("x"),  v_x,  kin->x );
        CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
        CheckBins( BinScheme("y"),  v_y,  kin->y );

        histSetFillList.clear();
        for(int bpt : v_pt) {
          for(int bz : v_z) {
            for(int bx : v_x) {
              for(int bq : v_q) {
                for(int by : v_y) {
                  histSetFillList.push_back(histSetJets[bpt][bz][bx][bq][by]);
                };
              };
            };
          };
        };

        for(Histos *H : histSetFillList) {
          H->Hist("pT_jet")->Fill(kin->pTjet,wJet);
          H->Hist("mT_jet")->Fill(jet.mt(),wJet);
          H->Hist("zJet")->Fill(kin->zjet,wJet);
          H->Hist("eta_jet")->Fill(jet.eta(),wJet);
          H->Hist("qT_jet")->Fill(kin->qTjet,wJet);
          H->Hist("qTQ")->Fill(kin->qTjet/sqrt(kin->Q2),wJet);
          for(int j = 0; j < kin->jperp.size(); j++) {
            H->Hist("jperp")->Fill(kin->jperp[j],wJet);
          };
        };

      };
    };

    #if INCCENTAURO == 1
    for(int brec = 0; brec < NrecMethodBins; brec++){
      TString recname = recMethodName.find(brec)->second;
      kin->CalculateDIS(recname);
      kin->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);

      if(kin->CutDIS()){
        for(int i = 0; i < kin->breitJetsRec.size(); i++){
          PseudoJet jet = kin->breitJetsRec[i];
          kin->CalculateBreitJetKinematics(jet);

          CheckBins( BinScheme("ptJet"), v_pt, kin->pTjet );
          CheckBins( BinScheme("zJet"), v_z, kin->zjet );
          CheckBins( BinScheme("x"),  v_x,  kin->x );
          CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
          CheckBins( BinScheme("y"),  v_y,  kin->y );

          histSetFillList.clear();
          for(int bpt : v_pt) {
            for(int bz : v_z) {
              for(int bx : v_x) {
                for(int bq : v_q) {
                  for(int by : v_y) {
                    histSetFillList.push_back(histSetBreitJets[bpt][bz][bx][bq][by][brec]);
                  };
                };
              };
            };
          };

          for(Histos *H : histSetFillList) {
            H->Hist("pT_jet")->Fill(kin->pTjet);
            H->Hist("mT_jet")->Fill(jet.mt());
            H->Hist("zJet")->Fill(kin->zjet);
            H->Hist("eta_jet")->Fill(jet.eta());
            H->Hist("qT_jet")->Fill(kin->qTjet);
            H->Hist("qTQ")->Fill(kin->qTjet/sqrt(kin->Q2));
            for(int j = 0; j < kin->jperp.size(); j++){
              H->Hist("jperp")->Fill(kin->jperp[j]);
            };
          };
        };
      };
    };
    #endif
    */


  };
  cout << "end event loop" << endl;
  // event loop end =========================================================

  // reset HD, to clean up after the event loop
  HD->ActivateAllNodes();
  HD->ClearOps();

  // calculate integrated luminosity
  Double_t lumi = wTotal/xsecTot; // [nb^-1]
  cout << "Integrated Luminosity:       " << lumi << "/nb" << endl;
  cout << sep << endl;

  // calculate cross sections, and print yields
  HD->Initial([&sep](){ cout << sep << endl << "Histogram Entries:" << endl; });
  HD->Final([&sep](){ cout << sep << endl; });
  HD->ForEach([&lumi](Histos *H){
    cout << H->GetSetTitle() << " ::: "
         << H->Hist("Q2vsX")->GetEntries()
         << endl;
    // calculate cross sections
    H->Hist("Q_xsec")->Scale(1./lumi); // TODO: generalize (`if (name contains "xsec") ...`)
  });
  HD->ExecuteAndClearOps();


  // write histograms
  cout << sep << endl;
  outFile->cd();
  if(writeSimpleTree) ST->WriteTree();
  HD->ForEach([this](Histos *H){ H->WriteHists(outFile); }); HD->ExecuteAndClearOps();
  HD->ForEach([this](Histos *H){ H->Write(); }); HD->ExecuteAndClearOps();
  /*
  for(Histos *H : histSetListJets) H->WriteHists(outFile);
  for(Histos *H : histSetListJets) H->Write();
  #if INCCENTAURO == 1
  for(Histos *H : histSetListBreitJets) H->WriteHists(outFile);
  for(Histos *H : histSetListBreitJets) H->Write();
  #endif
  */
  // write binning schemes
  for(auto const &kv : binSchemes) kv.second->Write(kv.first+"_bins"); // TODO: might not need this after DAG upgrade

  // close output
  outFile->Close();
  cout << outfileName << " written." << endl;

};



//=============================================



// access bin scheme by name
BinSet *Analysis::BinScheme(TString varname) {
  BinSet *ret;
  try { ret = binSchemes.at(varname); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: bin scheme "
         << varname << " not found" << endl;
    return nullptr;
  };
  return ret;
};

// add a new bin scheme
void Analysis::AddBinScheme(TString varname) {
  TString vartitle;
  try { vartitle = availableBinSchemes.at(varname); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: bin scheme "
         << varname << " not available... skipping..." << endl;
    return;
  };
  BinSet *B = new BinSet(varname,vartitle);
  binSchemes.insert(std::pair<TString,BinSet*>(varname,B));
};

// add a final state bin
void Analysis::AddFinalState(TString finalStateN, TString finalStateT, Int_t pid_) {
  // get bin number (we are adding a new bin, so new bin number = curent number of bins)
  Int_t binNum = BinScheme("finalState")->GetNumBins();
  // map : pid_ -> bin number
  PIDtoEnum.insert(std::pair<int,int>( pid_, binNum ));
  // map : bin number -> final state name (needed because this isn't stored in `CutDef`)
  finalStateName.insert(std::pair<int,TString>( binNum, finalStateN ));
  // build bin with custom `CutDef` ("custom" means that `CutDef` will not apply cuts,
  // rather the cuts are applied here)
  BinScheme("finalState")->BuildCustomBin(finalStateT);
};

// add reconstruction method bin
void Analysis::AddRecMethod(TString recMethodN, TString recMethodT){
  // creating binning for reconstruction methods, based on above function
  // for final states.
  Int_t binNum = BinScheme("recMethod")->GetNumBins();
  recMethodName.insert(std::pair<int,TString>(binNum, recMethodN));
  BinScheme("recMethod")->BuildCustomBin(recMethodT);
};


// return true, if a diagonal mode is on and this is an
// off-diagonal bin; if a diagonal mode is not on, always
// return false
Bool_t Analysis::CheckDiagonal(int cpt, int cx, int cz, int cq) {
  // TODO
  return true;
  /*
  if(diagonalPtXZ) return ( cpt!=cx || cx!=cz );
  else if(diagonalXZQ) return ( cx!=cz || cz!=cq );
  else return false;
  */
};

// scan through bin set `bs`, checking each one; the vector `v` will
// contain the list of array indices for which the cut on `var` is satisfied
/*
void Analysis::CheckBins(BinSet *bs, std::vector<int> &v, Double_t var) {
  v.clear();
  for(int b=0; b<bs->GetNumBins(); b++) {
    if(bs->Cut(b)->CheckCut(var)) v.push_back(b);
  };
};
*/

// get name of Histos object for specified bin
TString Analysis::GetHistosName(int cpt, int cx, int cz, int cq, int cy, int cfs) {
  TString retStr;
  retStr = "histos_";
  retStr += finalStateName[cfs];
  retStr += Form("_pt%d",cpt);
  retStr += Form("_x%d",cx);
  retStr += Form("_z%d",cz);
  retStr += Form("_q%d",cq);
  retStr += Form("_y%d",cy);
  return retStr;
};
TString Analysis::GetHistosTitle(int cpt, int cx, int cz, int cq, int cy, int cfs) {
  TString retStr;
  retStr  =        BinScheme("finalState")->Cut(cfs)->GetCutTitle();
  retStr += ", " + BinScheme("pt")->Cut(cpt)->GetCutTitle();
  retStr += ", " + BinScheme("x")->Cut(cx)->GetCutTitle();
  retStr += ", " + BinScheme("z")->Cut(cz)->GetCutTitle();
  retStr += ", " + BinScheme("q")->Cut(cq)->GetCutTitle();
  retStr += ", " + BinScheme("y")->Cut(cy)->GetCutTitle();
  return retStr;
};

TString Analysis::GetHistosNameJets(int cpt, int cz, int cx, int cq, int cy) {
  TString retStr;
  retStr = "histosJets_";
  retStr += Form("_pt_jet%d",cpt);
  retStr += Form("_z_jet%d",cz);
  retStr += Form("_x%d",cx);
  retStr += Form("_q%d",cq);
  retStr += Form("_y%d",cy);
  return retStr;
};

TString Analysis::GetHistosTitleJets(int cpt, int cz, int cx, int cq, int cy){
  TString retStr;
  retStr  =        BinScheme("ptJet")->Cut(cpt)->GetCutTitle();
  retStr += ", " + BinScheme("zJet")->Cut(cz)->GetCutTitle();
  retStr += ", " + BinScheme("x")->Cut(cx)->GetCutTitle();
  retStr += ", " + BinScheme("q")->Cut(cq)->GetCutTitle();
  retStr += ", " + BinScheme("y")->Cut(cy)->GetCutTitle();
  return retStr;
};
TString Analysis::GetHistosNameBreitJets(int cpt, int cz, int cx, int cq, int cy, int crec) {
  TString retStr;
  retStr = "histosBreitJets_";
  retStr += recMethodName[crec];
  retStr += Form("_pt_jet%d",cpt);
  retStr += Form("_z_jet%d",cz);
  retStr += Form("_x%d",cx);
  retStr += Form("_q%d",cq);
  retStr += Form("_y%d",cy);
  return retStr;
};

TString Analysis::GetHistosTitleBreitJets(int cpt, int cz, int cx, int cq, int cy, int crec){
  TString retStr;
  retStr  =        BinScheme("recMethod")->Cut(crec)->GetCutTitle();
  retStr +=  " " + BinScheme("ptJet")->Cut(cpt)->GetCutTitle();
  retStr += ", " + BinScheme("zJet")->Cut(cz)->GetCutTitle();
  retStr += ", " + BinScheme("x")->Cut(cx)->GetCutTitle();
  retStr += ", " + BinScheme("q")->Cut(cq)->GetCutTitle();
  retStr += ", " + BinScheme("y")->Cut(cy)->GetCutTitle();
  return retStr;
};

// destructor
Analysis::~Analysis() {
};

