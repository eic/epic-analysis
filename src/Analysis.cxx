#include "Analysis.h"

ClassImp(Analysis)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
Analysis::Analysis(
  TString infile_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_
)
  : infile(infile_)
  , eleBeamEn(eleBeamEn_)
  , ionBeamEn(ionBeamEn_)
  , crossingAngle(crossingAngle_)
{

  // set bin schemes
  AddBinScheme("pt","p_{T}");
  AddBinScheme("z","z");
  AddBinScheme("x","x");
  AddBinScheme("q","Q");
  AddBinScheme("y","y");
  AddBinScheme("pt_jet", "jet p_{T}");
  AddBinScheme("z_jet", "jet z");
  // final state bins (e.g., tracks or jets)
  AddBinScheme("finalState","finalState");
  AddFinalState("pipTrack","#pi^{+} tracks", 211);
  //AddFinalState("pimTrack","#pi^{-} tracks",-211);
  AddBinScheme("recMethod", "recMethod");
  AddRecMethod("Ele", "electron method");
  AddRecMethod("DA", "DA method");
  AddRecMethod("JB", "JB method");
  // initialize diagonalizer settings
  // TODO: generalized diagonalizer
  diagonalPtXZ = false;
  diagonalXZQ = false;
};


//=============================================
// perform the analysis
//=============================================
void Analysis::Execute() {

  // read delphes tree
  cout << "-- running analysis of " << infile << endl;
  TChain *chain = new TChain("Delphes");
  chain->Add(infile);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  Long64_t ENT = tr->GetEntries();

  // define output file
  TString outfileN = infile;
  outfileN(TRegexp("^.*/")) = "";
  outfileN = "out/histos."+outfileN;
  TFile *outfile = new TFile(outfileN,"RECREATE");

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itPIDSystemsTrack(tr->UseBranch("PIDSystemsTrack"));

  // number of bins
  const Int_t NptBins = BinScheme("pt")->GetNumBins();
  const Int_t NxBins = BinScheme("x")->GetNumBins();
  const Int_t NzBins = BinScheme("z")->GetNumBins();
  const Int_t NqBins = BinScheme("q")->GetNumBins();
  const Int_t NyBins = BinScheme("y")->GetNumBins();
  const Int_t NfinalStateBins = BinScheme("finalState")->GetNumBins();
  const Int_t NptjetBins = BinScheme("pt_jet")->GetNumBins();
  const Int_t NzjetBins = BinScheme("z_jet")->GetNumBins();
  const Int_t NrecMethodBins = BinScheme("recMethod")->GetNumBins();
  
  // sets of histogram sets
  // - `histSet` is a data structure for storing and organizing pointers to
  //   sets of histograms (`Histos` objects)
  // - `histSet*List` are used as temporary lists of relevant `Histos` pointers
  // - TODO: if we add one more dimension, 7D array will probably break; need
  //         better data structure
  Histos *histSet[NptBins][NxBins][NzBins][NqBins][NyBins][NfinalStateBins];
  Histos *histSetJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins];
  Histos *histSetBreitJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins][NrecMethodBins];
 
  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetListJets;
  std::vector<Histos*> histSetListBreitJets;

  std::vector<Histos*> histSetFillList;
  std::vector<int> v_pt, v_x, v_z, v_q, v_y;
  // instantiate Histos sets, and populate 
  TString histosN,histosT;

  cout << "Define track histograms..." << endl;
  for(int bpt=0; bpt<NptBins; bpt++) { // - loop over pT bins
    for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
      for(int bz=0; bz<NzBins; bz++) { // - loop over z bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q bins
          if(CheckDiagonal(bpt,bx,bz,bq)) continue; // diagonalizer
          for(int by=0; by<NyBins; by++) { // - loop over y bins
            for(int bfs=0; bfs<NfinalStateBins; bfs++) { // - loop over final states

              // set Histos name and title
              histosN = this->GetHistosName (bpt,bx,bz,bq,by,bfs);
              histosT = this->GetHistosTitle(bpt,bx,bz,bq,by,bfs);

              // define set of histograms for this bin
              histSet[bpt][bx][bz][bq][by][bfs] = new Histos(histosN,histosT);
              HS = histSet[bpt][bx][bz][bq][by][bfs]; // shorthand pointer

              // HISTOGRAMS ================================================
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

              // store cut definitions with histogram sets
              HS->AddCutDef(BinScheme("pt")->Cut(bpt));
              HS->AddCutDef(BinScheme("x")->Cut(bx));
              HS->AddCutDef(BinScheme("z")->Cut(bz));
              HS->AddCutDef(BinScheme("q")->Cut(bq));
              HS->AddCutDef(BinScheme("y")->Cut(by));
              HS->AddCutDef(BinScheme("finalState")->Cut(bfs));

              // add histogram set full list
              histSetList.push_back(histSet[bpt][bx][bz][bq][by][bfs]);
            };
          };
        };
      };
    };
  };

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
	    HS->DefineHist1D("z_jet","z","GeV", NBINS,0, 1);
	    HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
	    HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0); 
	    HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
	    HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);
	    // store cut definitions with histogram sets, then add histogram sets full list
	    HS->AddCutDef(BinScheme("pt_jet")->Cut(bpt));
	    HS->AddCutDef(BinScheme("z_jet")->Cut(bz));
	    HS->AddCutDef(BinScheme("x")->Cut(bx));	    
	    HS->AddCutDef(BinScheme("q")->Cut(bq));
	    HS->AddCutDef(BinScheme("y")->Cut(by));
	    histSetListJets.push_back(histSetJets[bpt][bz][bx][bq][by]);	  
	  };	
	};
      };
    };
  };
  #ifdef __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
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
	      HS->DefineHist1D("z_jet","z","GeV", NBINS,0, 1);
	      HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
	      HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0);
	      HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
	      HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);

	      // store cut definitions with histogram sets, then add histogram sets full list                                                                                                       
	      HS->AddCutDef(BinScheme("pt_jet")->Cut(bpt));
	      HS->AddCutDef(BinScheme("z_jet")->Cut(bz));
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
  // calculate integrated luminosity
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
  Double_t lumi = numGen/xsecTot; // [nb^-1]
  TString sep = "--------------------------------------------";
  cout << sep << endl;
  cout << "assumed total cross section: " << xsecTot << " nb" << endl;
  cout << "number of generated events:  " << numGen << endl;
  cout << "Integrated Luminosity:       " << lumi << "/nb" << endl;
  cout << sep << endl;


  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);


  // vars
  Double_t eleP,maxEleP;
  int pid,bFinalState;


  // event loop =========================================================
  ENT = 1000; // limiter
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

    // get hadronic final state variables
    kin->GetHadronicFinalState(itTrack, itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itPIDSystemsTrack, itParticle);
    // get vector of jets
    // should this have an option for clustering method?
    kin->CalculateDISbyElectron();
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
      kin->CalculateHadronKinematics();

      // apply cuts
      if(kin->CutFull()) {

        // decide which histogram sets to fill
        // -- `v_A` will be the list of bins that this event's variable `A`
        //    will be a part of; we use these lists to determine the list
        //    of histogram sets to fill
        // - check pT bin
        CheckBins( BinScheme("pt"), v_pt, kin->pT );
        CheckBins( BinScheme("x"),  v_x,  kin->x );
        CheckBins( BinScheme("z"),  v_z,  kin->z );
        CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
        CheckBins( BinScheme("y"),  v_y,  kin->y );


        // build list of histogram sets to fill
        histSetFillList.clear();
        for(int bpt : v_pt) {
          for(int bx : v_x) {
            for(int bz : v_z) {
              for(int bq : v_q) {
                for(int by : v_y) {
                  if(!CheckDiagonal(bpt,bx,bz,bq)) {
                    histSetFillList.push_back(histSet[bpt][bx][bz][bq][by][bFinalState]);
                  };
                };
              };
            };
          };
        };

        // loop through list of histogram sets, and fill them
        for(Histos *H : histSetFillList) {
          // DIS kinematics
          H->Hist("Q2vsX")->Fill(kin->x,kin->Q2);
          H->Hist("Q")->Fill(TMath::Sqrt(kin->Q2));
          H->Hist("x")->Fill(kin->x);
          H->Hist("W")->Fill(kin->W);
          H->Hist("y")->Fill(kin->y);
          // hadron 4-momentum
          H->Hist("pLab")->Fill(kin->pLab);
          H->Hist("pTlab")->Fill(kin->pTlab);
          H->Hist("etaLab")->Fill(kin->etaLab);
          H->Hist("phiLab")->Fill(kin->phiLab);
          // hadron kinematics
          H->Hist("z")->Fill(kin->z);
          H->Hist("pT")->Fill(kin->pT);
          H->Hist("qT")->Fill(kin->qT);
          H->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2));
          H->Hist("mX")->Fill(kin->mX);
          H->Hist("phiH")->Fill(kin->phiH);
          H->Hist("phiS")->Fill(kin->phiS);
          // cross sections
          H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),1.0/lumi);
        };

      };
    };

    // jet loop            
    if(kin->CutDIS()){
      for(int i = 0; i < kin->jetsRec.size(); i++){
	PseudoJet jet = kin->jetsRec[i];
	kin->CalculateJetKinematics(jet);
	
	// following same procedure as in track loop	
	CheckBins( BinScheme("pt_jet"), v_pt, kin->pTjet );
	CheckBins( BinScheme("z_jet"), v_z, kin->zjet );
	CheckBins( BinScheme("x"),  v_x,  kin->x );        
	CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
	CheckBins( BinScheme("y"),  v_y,  kin->y );
	
	histSetFillList.clear();
	for(int bpt : v_pt) {
	  for(int bz : v_z){
	    for(int bx : v_x) {      
	      for(int bq : v_q) {
		for(int by : v_y) {
		  
		  histSetFillList.push_back(histSetJets[bpt][bz][bx][bq][by]);
		  
		};
	      };
	      
	    };
	  };
	  for(Histos *H : histSetFillList) {	  
	    H->Hist("pT_jet")->Fill(kin->pTjet);
	    H->Hist("mT_jet")->Fill(jet.mt());
	    H->Hist("z_jet")->Fill(kin->zjet);
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
    #ifdef __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
    for(int brec = 0; brec < NrecMethodBins; brec++){
      TString recname = recMethodName.find(brec)->second;      
      kin->CalculateDIS(recname);      
      kin->GetBreitFrameJets(itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itParticle);
      
      if(kin->CutDIS()){
	for(int i = 0; i < kin->breitJetsRec.size(); i++){
	  PseudoJet jet = kin->breitJetsRec[i];
	  kin->CalculateBreitJetKinematics(jet);
	  
	  CheckBins( BinScheme("pt_jet"), v_pt, kin->pTjet );
	  CheckBins( BinScheme("z_jet"), v_z, kin->zjet );
	  CheckBins( BinScheme("x"),  v_x,  kin->x );
	  CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
	  CheckBins( BinScheme("y"),  v_y,  kin->y );
	  
	  histSetFillList.clear();
	  for(int bpt : v_pt) {
	    for(int bz : v_z){
	      for(int bx : v_x) {
		for(int bq : v_q) {
		  for(int by : v_y) {
		    
		    histSetFillList.push_back(histSetBreitJets[bpt][bz][bx][bq][by][brec]);
		    
		  };
		};
		
	      };
	    };
	    for(Histos *H : histSetFillList) {
	      H->Hist("pT_jet")->Fill(kin->pTjet);
	      H->Hist("mT_jet")->Fill(jet.mt());
	      H->Hist("z_jet")->Fill(kin->zjet);
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
      
    };
    #endif
  };
    cout << "end event loop" << endl;
  // event loop end =========================================================



  // print yields in each bin
  cout << sep << endl << "Histogram Entries:" << endl;
  for(Histos *H : histSetList) {
    cout << H->GetSetTitle() << " ::: "
         << H->Hist("Q2vsX")->GetEntries()
         << endl;
  };



  // write histograms
  cout << sep << endl;
  outfile->cd();
  for(Histos *H : histSetList) H->WriteHists(outfile);
  for(Histos *H : histSetList) H->Write();
  for(Histos *H : histSetListJets) H->WriteHists(outfile);
  for(Histos *H : histSetListJets) H->Write();
  #ifdef __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
  for(Histos *H : histSetListBreitJets) H->WriteHists(outfile);
  for(Histos *H : histSetListBreitJets) H->Write();
  #endif
  // write binning schemes
  for(auto const &kv : binSchemes) kv.second->Write(kv.first+"_bins");

  // close output
  outfile->Close();
  cout << outfileN << " written." << endl;

  // call draw program
  /*
  TString cmd = "./draw.exe "+outfileN;
  system(cmd.Data());
  */
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
void Analysis::AddBinScheme(TString varname, TString vartitle) {
  binSchemes.insert(
    std::pair<TString,BinSet*>(varname,new BinSet(varname,vartitle))
    );
  // TODO: for now, we need to have at least one bin in each dimension,
  // otherwise for loops won't run; when we generalize the `histSet` data
  // structure, hopefully we can also drop this requirement; the current
  // workaround is to add a `full` bin to each dimension
  if(varname!="finalState" && varname!="recMethod") BinScheme(varname)->BuildBin("Full");
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
  if(diagonalPtXZ) return ( cpt!=cx || cx!=cz );
  else if(diagonalXZQ) return ( cx!=cz || cz!=cq );
  else return false;
};

// scan through bin set `bs`, checking each one; the vector `v` will
// contain the list of array indices for which the cut on `var` is satisfied
void Analysis::CheckBins(BinSet *bs, std::vector<int> &v, Double_t var) {
  v.clear();
  for(int b=0; b<bs->GetNumBins(); b++) {
    if(bs->Cut(b)->CheckCut(var)) v.push_back(b);
  };
};

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
  retStr  =        BinScheme("pt_jet")->Cut(cpt)->GetCutTitle();
  retStr += ", " + BinScheme("z_jet")->Cut(cz)->GetCutTitle();
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
  retStr +=  " " + BinScheme("pt_jet")->Cut(cpt)->GetCutTitle();
  retStr += ", " + BinScheme("z_jet")->Cut(cz)->GetCutTitle();
  retStr += ", " + BinScheme("x")->Cut(cx)->GetCutTitle();
  retStr += ", " + BinScheme("q")->Cut(cq)->GetCutTitle();
  retStr += ", " + BinScheme("y")->Cut(cy)->GetCutTitle();
  return retStr;
};

// destructor
Analysis::~Analysis() {
};

