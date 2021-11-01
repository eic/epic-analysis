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
  entriesTot = 0;
};


// input files
//------------------------------------
// add a single file
bool Analysis::AddFile(TString fileName, Long64_t entries, Double_t xs, Double_t Q2min) {
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
  inXsecsTot.insert(inXsecsTot.begin() + insertIdx, xs);
  inQ2mins.insert(inQ2mins.begin() + insertIdx, Q2min);
  inEntries.insert(inEntries.begin() + insertIdx, entries);
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
    Long64_t entries;
    Double_t xs, Q2min;
    ss >> fileName >> entries >> xs >> Q2min;
    if (entries <= 0) {
      TFile* file = TFile::Open(fileName.c_str());
      if (file->IsZombie()) {
        cerr << "ERROR: Couldn't open input file '" << fileName << "'" << endl;
        return;
      }
      TTree* tree = file->Get<TTree>("Delphes");
      if (tree == nullptr) tree = file->Get<TTree>("events");
      if (tree == nullptr) {
        cerr << "ERROR: Couldn't find Delphes or events tree in file '" << fileName << "'" << endl;
        return;
      }
      entries = tree->GetEntries();
    }
    if (!ss) {
      continue;
    }
    if (!AddFile(TString(fileName.c_str()), entries, xs, Q2min)) {
      cerr << "ERROR: Couldn't add file '" << fileName << "'" << endl;
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
    HS->DefineHist1D("x_Res","x-x_{true}","", NBINS, -0.5, 0.5);
    HS->DefineHist1D("y_Res","y-y_{true}","", NBINS, -0.2, 0.2);
    HS->DefineHist1D("Q2_Res","Q2-Q2_{true}","GeV", NBINS, -0.5, 0.5);
    HS->DefineHist1D("phiH_Res","#phi_{h}-#phi_{h}^{true}","", NBINS, -TMath::Pi(), TMath::Pi());
    HS->DefineHist1D("phiS_Res","#phi_{S}-#phi_{S}^{true}","", NBINS, -0.1*TMath::Pi(), 0.1*TMath::Pi());
    HS->DefineHist1D("pT_Res","pT-pT^{true}","GeV", NBINS, -1.5, 1.5);
    HS->DefineHist2D("Q2vsXtrue","x","Q^{2}","","GeV^{2}",
        20,1e-4,1,
        10,1,1e4,
        true,true
        );
    HS->DefineHist2D("Q2vsXpurity","x","Q^{2}","","GeV^{2}",
        20,1e-4,1,
        10,1,1e4,
        true,true
        );
    HS->DefineHist2D("Q2vsX_zres","x","Q^{2}","","GeV^{2}",
        20,1e-4,1,
        10,1,1e4,
        true,true
        );
    HS->DefineHist2D("Q2vsX_pTres","x","Q^{2}","","GeV^{2}",
        20,1e-4,1,
        10,1,1e4,
        true,true
        );
    HS->DefineHist2D("Q2vsX_phiHres","x","Q^{2}","","GeV^{2}",
        20,1e-4,1,
        10,1,1e4,
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

void Analysis::CalculateEventQ2Weights() {
  Q2weights.resize(infiles.size());
  entriesTot = 0;
  for (Long64_t entry : inEntries) {
    entriesTot += entry;
  }
  xsecTot = inXsecsTot.front();
  for (Int_t idx = 0; idx < infiles.size(); ++idx) {
    Double_t xsecFactor = inXsecs[idx] / xsecTot;
    Double_t entries = 0.;
    for (Int_t idxC = 0; idxC <= idx; ++idxC) {
      // estimate how many events from this file lie in the given Q2 range
      entries += inEntries[idxC] * (inXsecs[idx] / inXsecsTot[idxC]);
    }
    Double_t numFactor = entries / entriesTot;
    Q2weights[idx] = xsecFactor / numFactor;
    cout << "Weights for " << infiles[idx] << ": " << numFactor << ", " << xsecFactor << endl;
  }
}

Double_t Analysis::GetEventQ2Weight(Double_t Q2, Int_t guess) {
  Int_t idx = GetEventQ2Idx(Q2, guess);
  if (idx == -1) {
    return 0.;
  } else {
    return Q2weights[idx];
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
  HD->Payload([&lumi](Histos *H){
    cout << H->GetSetTitle() << " ::: "
         << H->Hist("Q2vsX")->GetEntries()
         << endl;
    // calculate cross sections
    H->Hist("Q_xsec")->Scale(1./lumi); // TODO: generalize (`if (name contains "xsec") ...`)
    // divide resolution plots by true counts per x-Q2 bin
    H->Hist("Q2vsXpurity")->Divide(H->Hist("Q2vsXtrue"));
    H->Hist("Q2vsX_zres")->Divide(H->Hist("Q2vsXtrue"));
    H->Hist("Q2vsX_pTres")->Divide(H->Hist("Q2vsXtrue"));
    H->Hist("Q2vsX_phiHres")->Divide(H->Hist("Q2vsXtrue"));        
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
    H->Hist("pT_Res")->Fill( kin->pT - kinTrue->pT, wTrack );
    dynamic_cast<TH2*>(H->Hist("Q2vsXtrue"))->Fill(kinTrue->x,kinTrue->Q2,wTrack);
    if(kinTrue->z!=0) dynamic_cast<TH2*>(H->Hist("Q2vsX_zres"))->Fill(
      kinTrue->x,kinTrue->Q2,wTrack*( fabs(kinTrue->z - kin->z)/(kinTrue->z) ) );
    if(kinTrue->pT!=0) dynamic_cast<TH2*>(H->Hist("Q2vsX_pTres"))->Fill(
      kinTrue->x,kinTrue->Q2,wTrack*( fabs(kinTrue->pT - kin->pT)/(kinTrue->pT) ) );
    dynamic_cast<TH2*>(H->Hist("Q2vsX_phiHres"))->Fill(kinTrue->x,kinTrue->Q2,wTrack*( fabs(Kinematics::AdjAngle(kinTrue->phiH - kin->phiH) ) ) );
    
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

