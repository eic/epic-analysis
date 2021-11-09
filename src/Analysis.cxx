#include "Analysis.h"

#include <fstream>
#include <string>
#include <sstream>

ClassImp(Analysis)

using std::map;
using std::vector;
using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;

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
  // available variables for binning
  // - availableBinSchemes is a map from variable name to variable title
  // - try to avoid using underscores in the variable name (they are okay in the title);
  //   convention is camel case, starting lowercase 
  /* DIS */
  availableBinSchemes.insert(std::pair<TString,TString>("x","x"));
  availableBinSchemes.insert(std::pair<TString,TString>("q2","Q^{2}"));
  availableBinSchemes.insert(std::pair<TString,TString>("w","W"));
  availableBinSchemes.insert(std::pair<TString,TString>("y","y"));
  /* single hadron */
  availableBinSchemes.insert(std::pair<TString,TString>("p","p"));
  availableBinSchemes.insert(std::pair<TString,TString>("eta","#eta"));
  availableBinSchemes.insert(std::pair<TString,TString>("pt","p_{T}")); // transverse to q, in ion rest frame
  availableBinSchemes.insert(std::pair<TString,TString>("ptLab","p_{T}^{lab}")); // transverse to xy plane, in lab frame
  availableBinSchemes.insert(std::pair<TString,TString>("z","z"));
  availableBinSchemes.insert(std::pair<TString,TString>("qT","q_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("qTq","q_{T}/Q"));
  availableBinSchemes.insert(std::pair<TString,TString>("mX","M_{X}"));
  availableBinSchemes.insert(std::pair<TString,TString>("xF","x_{F}"));
  availableBinSchemes.insert(std::pair<TString,TString>("phiH","#phi_{h}"));
  availableBinSchemes.insert(std::pair<TString,TString>("phiS","#phi_{S}"));
  availableBinSchemes.insert(std::pair<TString,TString>("tSpin","spin"));
  /* jets */
  availableBinSchemes.insert(std::pair<TString,TString>("ptJet", "jet p_{T}"));
  availableBinSchemes.insert(std::pair<TString,TString>("zJet", "jet z"));


  // available final states
  // - specify which final states you want to include using `AddFinalState(TString name)`
  // - if you specify none, default final state(s) will be chosen for you
  availableBinSchemes.insert(std::pair<TString,TString>("finalState","finalState"));
  AddBinScheme("finalState");
  // - finalState name (ID) -> title
  finalStateToTitle.insert(std::pair<TString,TString>("pipTrack","#pi^{+} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("pimTrack","#pi^{-} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("KpTrack","K^{+} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("KmTrack","K^{-} tracks"));
  finalStateToTitle.insert(std::pair<TString,TString>("jet","jets"));
  // - PID -> finalState ID
  PIDtoFinalState.insert(std::pair<int, TString>( 211,"pipTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>(-211,"pimTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>( 321,"KpTrack"));
  PIDtoFinalState.insert(std::pair<int, TString>(-321,"KmTrack"));


  // kinematics reconstruction methods
  // - choose one of these methods using `SetReconMethod(TString name)`
  // - if you specify none, a default method will be chosen
  reconMethodToTitle.insert(std::pair<TString,TString>("Ele","Electron method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("DA","Double Angle method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("JB","Jacquet-Blondel method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("Mixed","Mixed method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("Sigma","Sigma method"));
  reconMethodToTitle.insert(std::pair<TString,TString>("eSigma","eSigma method"));

  // common settings defaults
  // - these settings can be set at the macro level
  writeSimpleTree = false;
  maxEvents = 0;
  useBreitJets = false;

  weight = new WeightsUniform();
  weightJet = new WeightsUniform();

  // miscellaneous
  infiles.clear();
  numGen = 0;
};


// input files
//------------------------------------
// add a single file
bool Analysis::AddFile(TString fileName, Double_t xs, Double_t Q2min) {
  if(fileName=="") return false;
  cout << "-- running analysis of " << fileName << endl;
  // insert in order of Q2min
  std::size_t insertIdx = 0;
  for (std::size_t idx = 0; idx < infiles.size(); ++idx) {
    if (infiles[idx] == fileName) {
      cerr << "File " << fileName << " appears twice in config file" << endl;
      return false;
    }
    if (Q2min == inQ2mins[idx]) {
      cerr << "Q2min " << Q2min << "appears twice in config file" << endl;
      return false;
    } else if (Q2min < inQ2mins[idx]) {
      break;
    } else {
      insertIdx += 1;
    }
  }
  infiles.insert(infiles.begin() + insertIdx, fileName);
  inXsecs.insert(inXsecs.begin() + insertIdx, xs);
  inQ2mins.insert(inQ2mins.begin() + insertIdx, Q2min);
  inEntries.insert(inEntries.begin() + insertIdx, 0);
  if (insertIdx == 0) {
    xsecTot = xs;
  }
  // adjust cross-sections to account for overlapping regions
  for (std::size_t idx = infiles.size() - 1; idx > insertIdx; --idx) {
    inXsecs[insertIdx] -= inXsecs[idx];
  }
  if (insertIdx > 0) {
    inXsecs[insertIdx - 1] -= inXsecs[insertIdx];
  }
  for (std::size_t idx = 0; idx < inXsecs.size(); ++idx) {
    if (inXsecs[idx] < 0.) {
      cerr << "Cross-sections must strictly decrease with stricter Q2min cuts" << endl;
      return false;
    }
  }
  return true;
}


// prepare for the analysis
//------------------------------------
void Analysis::Prepare() {
  ifstream fin(infileName);
  string line;
  while (std::getline(fin, line)) {
    stringstream ss(line);
    string fileName;
    Double_t xs, Q2min;
    ss >> fileName >> xs >> Q2min;
    if (!ss) {
      continue;
    }
    if (!AddFile(TString(fileName.c_str()), xs, Q2min)) {
      return;
    }
  }
  if (infiles.empty()) {
    cerr << "ERROR: no input files have been specified" << endl;
    return;
  }

  // set output file name
  outfileName = "out/"+outfilePrefix+".root";

  // open output file
  cout << "-- output file: " << outfileName << endl;
  outFile = new TFile(outfileName,"RECREATE");

  // instantiate shared objects
  kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);
  kinTrue = new Kinematics(eleBeamEn, ionBeamEn, crossingAngle);
  ST = new SimpleTree("tree",kin,kinTrue);


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
  HD = new HistosDAG();
  HD->Build(binSchemes);


  // DEFINE HISTOGRAMS ------------------------------------
  HD->Payload([this](Histos *HS){
    // -- Full phase space histogram
    HS->DefineHist4D(
        "full_xsec",
        "x","Q^{2}","z","p_{T}",
        "","GeV^{2}","","GeV",
        NBINS_FULL,1e-3,1,
        NBINS_FULL,1,100,
        NBINS_FULL,0,1,
        NBINS_FULL,0,2,
        true,true
        );
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
    HS->DefineHist2D("phiHvsPhiS","#phi_{S}","#phi_{h}","","",
        25,-TMath::Pi(),TMath::Pi(),
        25,-TMath::Pi(),TMath::Pi());
    HS->DefineHist1D("phiSivers","#phi_{Sivers}","",NBINS,-TMath::Pi(),TMath::Pi());
    HS->DefineHist1D("phiCollins","#phi_{Collins}","",NBINS,-TMath::Pi(),TMath::Pi());
    HS->DefineHist2D("etaVsP","p","#eta","GeV","",
        NBINS,0.1,100,
        NBINS,-4,4,
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
    // -- resolutions
    HS->DefineHist1D("Q2_Res","Q2_{true}-Q2","GeV^{2}", NBINS, -2, 2);
    HS->DefineHist1D("x_Res","x_{true}-x","", NBINS, -2, 2);
    HS->DefineHist1D("y_Res","y_{true}-y","", NBINS, -2, 2);
    HS->DefineHist1D("z_Res","z_{true}-z","", NBINS, -2, 2);
    HS->DefineHist1D("pT_Res","pT_{true}-pT","GeV", NBINS, -2, 2);
    HS->DefineHist1D("phiH_Res","#phi_{h}^{true}-#phi_{h}","", NBINS, -TMath::Pi(), TMath::Pi());
    HS->DefineHist1D("phiS_Res","#phi_{S}^{true}-#phi_{S}","", NBINS, -TMath::Pi(), TMath::Pi());

    // resolutions vs. z on x axis.
    HS->DefineHist2D("z_Q2_Res","z","","#sigma_{Q2}","", NBINS, 0, 1, NBINSRES, -0.5, 0.5);//TODO: Fill these
    HS->DefineHist2D("z_x_Res","z","","#sigma_{x}","", NBINS, 0, 1, NBINSRES, -0.5, 0.5);
    HS->DefineHist2D("z_y_Res","z","","#sigma_{y}","", NBINS, 0, 1, NBINSRES, -0.5, 0.5);
    HS->DefineHist2D("z_z_Res","z","","#sigma_{z}","", NBINS, 0, 1, NBINSRES, -0.5, 0.5);
    HS->DefineHist2D("z_pT_Res","z","","#sigma_{pT}","", NBINS, 0, 1, NBINSRES, -0.5, 0.5);
    HS->DefineHist2D("z_phiH_Res","z","","#sigma_{#phi_{h}^{true}}","", NBINS, 0, 1, NBINSRES, -TMath::Pi(), TMath::Pi());
    HS->DefineHist2D("z_phiS_Res","z","","#sigma_{#phi_{S}^{true}}","", NBINS, 0, 1, NBINSRES, -TMath::Pi(), TMath::Pi());

    // 1D z-binned counts and summed resolutions
    HS->DefineHist1D("z_true","z","", NBINS, 0, 1);
    HS->DefineHist1D("z_purity","purity","", NBINS, 0, 1);

    // 2D Q2 vs. x binned resolutions
    HS->DefineHist2D("Q2vsXtrue","x","Q^{2}","","GeV^{2}",
        // 20,1e-4,1,//TODO: OLD -> Might revert...
        // 10,1,1e4,
        // true,true
        NBINS,1e-3,1,
        NBINS,1,100,
        true,true
        );


    // -- reconstructed vs. generated
    HS->DefineHist2D("x_RvG","generated x","reconstructed x","","",
        NBINS,1e-3,1,
        NBINS,1e-3,1,
        true,true
        );
    HS->DefineHist2D("phiH_RvG","generated #phi_{h}","reconstructed #phi_{h}","","",
        NBINS,-TMath::Pi(),TMath::Pi(),
        NBINS,-TMath::Pi(),TMath::Pi()
        );
    HS->DefineHist2D("phiS_RvG","generated #phi_{S}","reconstructed #phi_{S}","","",
        NBINS,-TMath::Pi(),TMath::Pi(),
        NBINS,-TMath::Pi(),TMath::Pi()
        );
  });
  HD->ExecuteAndClearOps();

  // initialize total weights
  wTrackTotal = 0.;
  wJetTotal = 0.;
};

void Analysis::CountEvent(Double_t Q2, Int_t guess) {
  Int_t idx = GetEventQ2Idx(Q2, guess);
  if (idx == -1) {
  } else {
    inEntries[idx] += 1;
    numGen += 1;
  }
}

Double_t Analysis::GetEventQ2Weight(Double_t Q2, Int_t guess) {
  Int_t idx = GetEventQ2Idx(Q2, guess);
  if (idx == -1) {
    return 0.;
  } else {
    Double_t xsecFactor = inXsecs[idx] / xsecTot;
    Double_t numFactor = (Double_t) inEntries[idx] / (Double_t) numGen;
    return xsecFactor / numFactor;
  }
}

Int_t Analysis::GetEventQ2Idx(Double_t Q2, Int_t guess) {
  Int_t idx = guess;
  if (Q2 < inQ2mins[idx]) {
    do {
      idx -= 1;
    } while (idx >= 0 && Q2 < inQ2mins[idx]);
    return idx;
  } else {
    while (idx + 1 < inQ2mins.size() && Q2 >= inQ2mins[idx + 1]) {
      idx += 1;
    }
    return idx;
  }
}

// finish the analysis
//-----------------------------------
void Analysis::Finish() {

  // reset HD, to clean up after the event loop
  HD->ActivateAllNodes();
  HD->ClearOps();

  // calculate integrated luminosity
  Double_t lumi = wTrackTotal/xsecTot; // [nb^-1]
  cout << "Integrated Luminosity:       " << lumi << "/nb" << endl;
  cout << sep << endl;

  // calculate cross sections, and print yields
  HD->Initial([this](){ cout << sep << endl << "Histogram Entries:" << endl; });
  HD->Final([this](){ cout << sep << endl; });
  HD->Payload([&lumi,this](Histos *H){
    cout << H->GetSetTitle() << " ::: "
         << H->Hist("Q2vsX")->GetEntries()
         << endl;
    // calculate cross sections
    H->Hist("Q_xsec")->Scale(1./lumi); // TODO: generalize (`if (name contains "xsec") ...`)

    // // Normalize and sqrt 1D z-binned resolution StdDevs
    H->Hist("z_purity")->Divide(H->Hist("z_true"));

  });

  HD->ExecuteAndClearOps();

  // write histograms
  cout << sep << endl;
  cout << "writing ROOT file..." << endl;
  outFile->cd();
  if(writeSimpleTree) ST->WriteTree();
  HD->Payload([this](Histos *H){ H->WriteHists(outFile); }); HD->ExecuteAndClearOps();
  HD->Payload([this](Histos *H){ H->Write(); }); HD->ExecuteAndClearOps();

  // write binning schemes
  for(auto const &kv : binSchemes) {
    if(kv.second->GetNumBins()>0) kv.second->Write("binset__"+kv.first);
  };

  // close output
  outFile->Close();
  cout << outfileName << " written." << endl;
};


// access bin scheme by name
//------------------------------------
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
//------------------------------------
void Analysis::AddBinScheme(TString varname) {
  TString vartitle;
  // TODO [low priority]: would be nice to make this lookup case insensitive
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
//------------------------------------
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


// access HistosDAG
//------------------------------------
HistosDAG *Analysis::GetHistosDAG() { return HD; };


// lambda to check which bins an observable is in, during DAG breadth
// traversal; it requires `finalStateID`, `valueMap`, and will
// activate/deactivate bin nodes accoding to values in `valuMap`
//--------------------------------------------------------------------
std::function<void(Node*)> Analysis::CheckBin() {
  return [this](Node *N){
    if(N->GetNodeType()==NT::bin) {
      Bool_t active;
      Double_t val;
      if(N->GetVarName()=="finalState") active = (N->GetCut()->GetCutID()==finalStateID);
      else {
        try {
          // get value associated to this variable, and check cut
          val = valueMap.at(N->GetVarName());
          active = N->GetCut()->CheckCut(val);
        } catch(const std::out_of_range &ex) {
          /* if this variable is not found in `valueMap`, then just activate
           * the node; this can happen if you are looking at jets AND tracks
           * final states, and you defined a binning scheme only valid for
           * tracks, but not for jets, e.g., `phiS`; if the current finalState
           * you are checking is a jet, we don't need to check phiS, so just
           * activate the node and ignore that cut
           */
          active = true;
        };
      };
      if(active) activeEvent=true;
      N->SetActiveState(active);
    };
  };
};


// FillHistos methods: check bins and fill associated histograms
// - checks which bins the track/jet/etc. falls in
// - fills the histograms in the associated Histos objects
//--------------------------------------------------------------------------
// tracks (single particles)
void Analysis::FillHistosTracks() {

  // add kinematic values to `valueMap`
  valueMap.clear();
  activeEvent = false;
  /* DIS */
  valueMap.insert(std::pair<TString,Double_t>( "x", kin->x ));
  valueMap.insert(std::pair<TString,Double_t>( "q2", kin->Q2 ));
  valueMap.insert(std::pair<TString,Double_t>( "w", kin->W ));
  valueMap.insert(std::pair<TString,Double_t>( "y", kin->y ));
  /* single hadron */
  valueMap.insert(std::pair<TString,Double_t>( "p", kin->pLab ));
  valueMap.insert(std::pair<TString,Double_t>( "eta", kin->etaLab ));
  valueMap.insert(std::pair<TString,Double_t>( "pt", kin->pT ));
  valueMap.insert(std::pair<TString,Double_t>( "ptLab", kin->pTlab ));
  valueMap.insert(std::pair<TString,Double_t>( "z", kin->z ));
  valueMap.insert(std::pair<TString,Double_t>( "qT", kin->qT ));
  valueMap.insert(std::pair<TString,Double_t>( "qTq", kin->qT/TMath::Sqrt(kin->Q2) ));
  valueMap.insert(std::pair<TString,Double_t>( "mX", kin->mX ));
  valueMap.insert(std::pair<TString,Double_t>( "xF", kin->xF ));
  valueMap.insert(std::pair<TString,Double_t>( "phiH", kin->phiH ));
  valueMap.insert(std::pair<TString,Double_t>( "phiS", kin->phiS ));
  valueMap.insert(std::pair<TString,Double_t>( "tSpin", (Double_t)kin->tSpin ));

  // check bins
  // - activates HistosDAG bin nodes which contain this track
  // - sets `activeEvent` if there is at least one multidimensional bin to fill
  HD->TraverseBreadth(CheckBin());
  if(!activeEvent) return;
  
  // fill histograms, for activated bins only
  HD->Payload([this](Histos *H){
    // Full phase space.
    H->Hist4("full_xsec")->Fill(kin->x,kin->Q2,kin->pT,kin->z,wTrack);
    // DIS kinematics
    dynamic_cast<TH2*>(H->Hist("Q2vsX"))->Fill(kin->x,kin->Q2,wTrack);
    H->Hist("Q")->Fill(TMath::Sqrt(kin->Q2),wTrack);
    H->Hist("x")->Fill(kin->x,wTrack);
    H->Hist("W")->Fill(kin->W,wTrack);
    H->Hist("y")->Fill(kin->y,wTrack);
    // hadron 4-momentum
    H->Hist("pLab")->Fill(kin->pLab,wTrack);
    H->Hist("pTlab")->Fill(kin->pTlab,wTrack);
    H->Hist("etaLab")->Fill(kin->etaLab,wTrack);
    H->Hist("phiLab")->Fill(kin->phiLab,wTrack);
    // hadron kinematics
    H->Hist("z")->Fill(kin->z,wTrack);
    H->Hist("pT")->Fill(kin->pT,wTrack);
    H->Hist("qT")->Fill(kin->qT,wTrack);
    if(kin->Q2!=0) H->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2),wTrack);
    H->Hist("mX")->Fill(kin->mX,wTrack);
    H->Hist("phiH")->Fill(kin->phiH,wTrack);
    H->Hist("phiS")->Fill(kin->phiS,wTrack);
    dynamic_cast<TH2*>(H->Hist("phiHvsPhiS"))->Fill(kin->phiS,kin->phiH,wTrack);
    H->Hist("phiSivers")->Fill(Kinematics::AdjAngle(kin->phiH - kin->phiS),wTrack);
    H->Hist("phiCollins")->Fill(Kinematics::AdjAngle(kin->phiH + kin->phiS),wTrack);
    dynamic_cast<TH2*>(H->Hist("etaVsP"))->Fill(kin->pLab,kin->etaLab,wTrack); // TODO: lab-frame p, or some other frame?
    // cross sections (divide by lumi after all events processed)
    H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),wTrack);
    // resolutions
    H->Hist("x_Res")->Fill( kin->x - kinTrue->x, wTrack );
    H->Hist("y_Res")->Fill( kin->y - kinTrue->y, wTrack );
    H->Hist("Q2_Res")->Fill( kin->Q2 - kinTrue->Q2, wTrack );
    H->Hist("phiH_Res")->Fill( Kinematics::AdjAngle(kin->phiH - kinTrue->phiH), wTrack );
    H->Hist("phiS_Res")->Fill( Kinematics::AdjAngle(kin->phiS - kinTrue->phiS), wTrack );

    // z binned resolutions
    if(kinTrue->Q2!=0) dynamic_cast<TH2*>(H->Hist("z_Q2_Res"))->Fill( kinTrue->z, (kinTrue->Q2 - kin->Q2)/kinTrue->Q2, wTrack );
    if(kinTrue->z!=0)  dynamic_cast<TH2*>(H->Hist("z_x_Res"))->Fill( kinTrue->z, (kinTrue->x - kin->x)/kinTrue->x, wTrack );
    if(kinTrue->y!=0)  dynamic_cast<TH2*>(H->Hist("z_y_Res"))->Fill( kinTrue->z, (kinTrue->y - kin->y)/kinTrue->y, wTrack );
    if(kinTrue->z!=0)  dynamic_cast<TH2*>(H->Hist("z_z_Res"))->Fill( kinTrue->z, (kinTrue->z - kin->z)/kinTrue->z, wTrack );
    if(kinTrue->pT!=0) dynamic_cast<TH2*>(H->Hist("z_pT_Res"))->Fill( kinTrue->z, (kinTrue->pT - kin->pT)/kinTrue->pT, wTrack );
    dynamic_cast<TH2*>(H->Hist("z_phiH_Res"))->Fill( kinTrue->z, Kinematics::AdjAngle(kin->phiH - kinTrue->phiH), wTrack );
    dynamic_cast<TH2*>(H->Hist("z_phiS_Res"))->Fill( kinTrue->z, Kinematics::AdjAngle(kin->phiS - kinTrue->phiS), wTrack );

    // purities
    H->Hist("z_true")->Fill(kinTrue->z, wTrack );
    if( (H->Hist("z_true"))->FindBin(kinTrue->z) == (H->Hist("z_true"))->FindBin(kin->z) ) H->Hist("z_purity")->Fill(kin->z,wTrack);

    if( (H->Hist("Q2vsXtrue"))->FindBin(kinTrue->x,kinTrue->Q2) == (H->Hist("Q2vsXtrue"))->FindBin(kin->x,kin->Q2) ) dynamic_cast<TH2*>(H->Hist("Q2vsXpurity"))->Fill(kin->x,kin->Q2,wTrack);
    
    // -- reconstructed vs. generated
    dynamic_cast<TH2*>(H->Hist("x_RvG"))->Fill(kinTrue->x,kin->x,wTrack);
    dynamic_cast<TH2*>(H->Hist("phiH_RvG"))->Fill(kinTrue->phiH,kin->phiH,wTrack);
    dynamic_cast<TH2*>(H->Hist("phiS_RvG"))->Fill(kinTrue->phiS,kin->phiS,wTrack);
  });
  // execute the payload
  // - save time and don't call `ClearOps` (next loop will overwrite lambda)
  // - called with `activeNodesOnly==true` since we only want to fill bins associated
  //   with this track
  HD->ExecuteOps(true);
};

// jets
void Analysis::FillHistosJets() {

  // add kinematic values to `valueMap`
  valueMap.clear();
  activeEvent = false;
  /* DIS */
  valueMap.insert(std::pair<TString,Double_t>(  "x",      kin->x      ));
  valueMap.insert(std::pair<TString,Double_t>(  "q2",     kin->Q2     ));
  valueMap.insert(std::pair<TString,Double_t>(  "y",      kin->y      ));
  /* jets */
  valueMap.insert(std::pair<TString,Double_t>(  "ptJet",  kin->pTjet  ));
  valueMap.insert(std::pair<TString,Double_t>(  "zJet",   kin->zjet   ));

  // check bins
  // - activates HistosDAG bin nodes which contain this track
  // - sets `activeEvent` if there is at least one multidimensional bin to fill
  HD->TraverseBreadth(CheckBin());
  if(!activeEvent) return;

  // fill histograms, for activated bins only
  HD->Payload([this](Histos *H){
    dynamic_cast<TH2*>(H->Hist("Q2vsX"))->Fill(kin->x,kin->Q2,wJet);
    // jet kinematics
    H->Hist("pT_jet")->Fill(kin->pTjet,wJet);
    H->Hist("mT_jet")->Fill(jet.mt(),wJet);
    H->Hist("z_jet")->Fill(kin->zjet,wJet);
    H->Hist("eta_jet")->Fill(jet.eta(),wJet);
    H->Hist("qT_jet")->Fill(kin->qTjet,wJet);
    if(kin->Q2!=0) H->Hist("qTQ_jet")->Fill(kin->qTjet/sqrt(kin->Q2),wJet);
    for(int j = 0; j < kin->jperp.size(); j++) {
      H->Hist("jperp")->Fill(kin->jperp[j],wJet);
    };
  });
  // execute the payload
  // - save time and don't call `ClearOps` (next loop will overwrite lambda)
  // - called with `activeNodesOnly==true` since we only want to fill bins associated
  //   with this jet
  HD->ExecuteOps(true);
};

// destructor
Analysis::~Analysis() {
};

