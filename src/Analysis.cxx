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
  , reconMethod("")
{
  // build list of variables available for binning (paired with titles)
  // - availableBinSchemes is a map from variable name to variable title
  // - try to avoid using underscores in the variable name (they are okay in the title)
  availableBinSchemes.insert(std::pair<TString,TString>("pt","p_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("p","p"));
  availableBinSchemes.insert(std::pair<TString,TString>("z","z"));
  availableBinSchemes.insert(std::pair<TString,TString>("x","x"));
  availableBinSchemes.insert(std::pair<TString,TString>("q2","Q^{2}"));
  availableBinSchemes.insert(std::pair<TString,TString>("y","y"));
  availableBinSchemes.insert(std::pair<TString,TString>("eta","#eta"));
  availableBinSchemes.insert(std::pair<TString,TString>("ptJet", "jet p_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("zJet", "jet z"));
  availableBinSchemes.insert(std::pair<TString,TString>("finalState","finalState"));
  availableBinSchemes.insert(std::pair<TString,TString>("recMethod","recMethod"));

  // available final states
  AddBinScheme("finalState");
  // - finalState name (ID) -> title
  finalStateToTitle.insert(std::pair<TString,TString>("pipTrack","#pi^{+} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("pimTrack","#pi^{-} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("KpTrack","K^{+} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("KmTrack","K^{-} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("jet","jets"));
  finalStateToTitle.insert(std::pair<TString,TString>("jetBreit","Breit jets"));
  // - PID -> finalState ID
  PIDtoFinalState.insert(std::pair<int, TString>( 211,"pipTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>(-211,"pimTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>( 321,"KpTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>(-321,"KmTrack"));

  // kinematics reconstruction methods
  reconMethodToTitle.insert(std::pair<TString,TString>("Ele","Electron method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("DA","Double Angle method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("JB","Jacquet-Blondel method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("Mixed","Mixed method"));

  // initializations
  writeSimpleTree = false;
  maxEvents = 0;
  useBreitJets = false;
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
  outfileName(TRegexp("\\.\\.")) = "."; // remove double dot
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

  // if there are no final states defined, default to definitions here:
  if(BinScheme("finalState")->GetNumBins()==0) {
    std::cout << "NOTE: adding pi+ tracks for final state, since you specified none" << std::endl;
    AddFinalState("pipTrack");
  };

  // if no reconstruction method is set, choose a default here
  if(reconMethod=="") {
    std::cout << "NOTE: no recon method specified, default to electron method" << std::endl;
    SetReconMethod("Ele");
  };

  // build HistosDAG with specified binning
  HD->Build(binSchemes);


  // DEFINE HISTOGRAMS ================================================
  HD->Payload([this](Histos *HS){
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
    HS->DefineHist2D("etaVsP","p","#eta","GeV","",
        NBINS,0.1,100,
        NBINS,-5,5,
        true,false
        );
    // -- single-hadron cross sections
    //HS->DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
    HS->DefineHist1D("Q_xsec","Q","GeV",10,1.0,10.0,true,true); // log
    HS->Hist("Q_xsec")->SetMinimum(1e-10);
    // -- jet kinematics
    HS->DefineHist1D("pT_jet","jet p_{T}","GeV", NBINS, 1e-2, 50);
    HS->DefineHist1D("mT_jet","jet m_{T}","GeV", NBINS, 1e-2, 20);
    HS->DefineHist1D("z_jet","jet z","", NBINS,0, 1);
    HS->DefineHist1D("eta_jet","jet #eta_{lab}","", NBINS,-5,5);
    HS->DefineHist1D("qT_jet","jet q_{T}", "GeV", NBINS, 0, 10.0);
    HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
    HS->DefineHist1D("qTQ_jet","jet q_{T}/Q","", NBINS, 0, 3.0);
  });
  HD->ExecuteAndClearOps();


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


  // lambda to check which bins an observable is in; it requires
  // `finalStateID`, `valueMap`, and will activate/deactivate bin nodes
  // accoding to values in `valuMap`
  auto checkBins = [this](Node *N){
    if(N->GetNodeType()==NT::bin) {
      Bool_t active;
      if(N->GetVarName()=="finalState") active = (N->GetCut()->GetCutID()==finalStateID);
      else {
        auto val = valueMap.at(N->GetVarName());
        active = N->GetCut()->CheckCut(val);
      };
      if(active) activeEvent=true;
      N->SetActiveState(active);
    };
  };


  // event loop =========================================================
  if(maxEvents>0) ENT = maxEvents; // limiter
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
    kin->CalculateDIS(reconMethod); // reconstructed
    kinTrue->CalculateDIS(reconMethod); // generated (truth)

    // get vector of jets
    // TODO: should this have an option for clustering method?
    kin->GetJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);

    // track loop
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
      
      // asymmetry injection
      //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

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
        activeEvent = false;
        HD->TraverseBreadth(checkBins);

        // fill histograms
        HD->Payload([this,&w](Histos *H){
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
          dynamic_cast<TH2*>(H->Hist("etaVsP"))->Fill(kin->pLab,kin->etaLab,w); // TODO: lab-frame p, or some other frame?
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
    }; // end track loop

    // jet loop
    finalStateID = "jet";
    if(activeFinalStates.find(finalStateID)!=activeFinalStates.end()) {

      #if INCCENTAURO == 1
      if(useBreitJets) kin->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
      #endif

      if(kin->CutDIS()){

        Double_t wJet = weightJet->GetWeight(*kin); // TODO: should we separate weights for breit and non-breit jets?
        wJetTotal += wJet;

        Int_t nJets;
        if(useBreitJets) nJets = kin->breitJetsRec.size();
        else      nJets = kin->jetsRec.size();

        for(int i = 0; i < kin->jetsRec.size(); i++){

          PseudoJet jet;
          if(useBreitJets) {
            #if INCCENTAURO == 1
            jet = kin->breitJetsRec[i];
            kin->CalculateBreitJetKinematics(jet);
            #endif
          } else {
            jet = kin->jetsRec[i];
            kin->CalculateJetKinematics(jet);
          };

          // map varNames to values
          // following same procedure as in track loop	
          valueMap.clear();
          valueMap.insert(std::pair<TString,Double_t>(  "ptJet",  kin->pTjet  ));
          valueMap.insert(std::pair<TString,Double_t>(  "zJet",   kin->zjet   ));
          valueMap.insert(std::pair<TString,Double_t>(  "x",      kin->x      ));
          valueMap.insert(std::pair<TString,Double_t>(  "q2",     kin->Q2     ));
          valueMap.insert(std::pair<TString,Double_t>(  "y",      kin->y      ));

          // check which bins the event falls in
          activeEvent = false;
          HD->TraverseBreadth(checkBins);

          // fill histograms
          HD->Payload([this,&wJet,&jet](Histos *H){
            // jet kinematics
            H->Hist("pT_jet")->Fill(kin->pTjet,wJet);
            H->Hist("mT_jet")->Fill(jet.mt(),wJet);
            H->Hist("z_jet")->Fill(kin->zjet,wJet);
            H->Hist("eta_jet")->Fill(jet.eta(),wJet);
            H->Hist("qT_jet")->Fill(kin->qTjet,wJet);
            H->Hist("qTQ_jet")->Fill(kin->qTjet/sqrt(kin->Q2),wJet);
            for(int j = 0; j < kin->jperp.size(); j++) {
              H->Hist("jperp")->Fill(kin->jperp[j],wJet);
            };
          });
          HD->ExecuteOps(true); // save time and don't ClearOps (next loop will overwrite lambda)

        };
      };
    }; // end jet loop

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
  HD->Payload([&lumi](Histos *H){
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
  HD->Payload([this](Histos *H){ H->WriteHists(outFile); }); HD->ExecuteAndClearOps();
  HD->Payload([this](Histos *H){ H->Write(); }); HD->ExecuteAndClearOps();

  // write binning schemes
  for(auto const &kv : binSchemes) kv.second->Write("binset__"+kv.first);

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
  if(binSchemes.find(varname)==binSchemes.end()) { // (duplicate prevention)
    BinSet *B = new BinSet(varname,vartitle);
    binSchemes.insert(std::pair<TString,BinSet*>(varname,B));
  };
};


// add a final state bin
void Analysis::AddFinalState(TString finalStateN) {
  TString finalStateT;
  try { finalStateT = finalStateToTitle.at(finalStateN); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: final state "
         << finalStateN << " not available... skipping..." << endl;
    return;
  };
  BinScheme("finalState")->BuildExternalBin(finalStateN,finalStateT);
  activeFinalStates.insert(finalStateN);
};


// destructor
Analysis::~Analysis() {
};

