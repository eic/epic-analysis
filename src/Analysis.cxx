#include "Analysis.h"

ClassImp(Analysis)

using std::cout;
using std::cerr;
using std::endl;

// constructor
Analysis::Analysis(
  TString infileName_,
  TString outfilePrefix_
)
  : infileName(infileName_)
  , outfilePrefix(outfilePrefix_)
  , reconMethod("")
  , finalStateID("")
  , eleBeamEn(10.0)
  , ionBeamEn(100.0)
  , crossingAngle(-25.0)
  , totalCrossSection(1e6)
{
  // available variables for binning
  // - availableBinSchemes is a map from variable name to variable title
  // - try to avoid using underscores in the variable name (they are okay in the title);
  //   convention is camel case, starting lowercase 
  /* DIS */
  availableBinSchemes.insert({ "x",     "x"           });
  availableBinSchemes.insert({ "q2",    "Q^{2}"       });
  availableBinSchemes.insert({ "w",     "W"           });
  availableBinSchemes.insert({ "y",     "y"           });
  /* single hadron */
  availableBinSchemes.insert({ "p",     "p"           });
  availableBinSchemes.insert({ "eta",   "#eta"        });
  availableBinSchemes.insert({ "pt",    "p_{T}"       }); // transverse to q, in ion rest frame
  availableBinSchemes.insert({ "ptLab", "p_{T}^{lab}" }); // transverse to xy plane, in lab frame
  availableBinSchemes.insert({ "z",     "z"           });
  availableBinSchemes.insert({ "qT",    "q_{T}"       });
  availableBinSchemes.insert({ "qTq",   "q_{T}/Q"     });
  availableBinSchemes.insert({ "mX",    "M_{X}"       });
  availableBinSchemes.insert({ "xF",    "x_{F}"       });
  availableBinSchemes.insert({ "phiH",  "#phi_{h}"    });
  availableBinSchemes.insert({ "phiS",  "#phi_{S}"    });
  availableBinSchemes.insert({ "tSpin", "spin"        });
  availableBinSchemes.insert({ "lSpin", "spinL"       });


  // available final states
  // - specify which final states you want to include using `AddFinalState(TString name)`
  // - if you specify none, default final state(s) will be chosen for you
  availableBinSchemes.insert({"finalState","finalState"});
  AddBinScheme("finalState");
  // - finalState name (ID) -> title
  finalStateToTitle.insert({ "pipTrack", "#pi^{+} tracks" });
  finalStateToTitle.insert({ "pimTrack", "#pi^{-} tracks" });
  finalStateToTitle.insert({ "KpTrack",  "K^{+} tracks"   });
  finalStateToTitle.insert({ "KmTrack",  "K^{-} tracks"   });
  finalStateToTitle.insert({ "pTrack",   "p^{+} tracks"   });
  // - PID -> finalState ID
  PIDtoFinalState.insert({  211,  "pipTrack" });
  PIDtoFinalState.insert({ -211,  "pimTrack" });
  PIDtoFinalState.insert({  321,  "KpTrack"  });
  PIDtoFinalState.insert({ -321,  "KmTrack"  });
  PIDtoFinalState.insert({  2212, "pTrack"   });

  // jets
#ifndef EXCLUDE_DELPHES
  // available variables for binning
  availableBinSchemes.insert({ "ptJet", "jet p_{T}" });
  availableBinSchemes.insert({ "zJet",  "jet z"     });
  // available final states
  finalStateToTitle.insert({ "jet", "jets" });
#endif

  // kinematics reconstruction methods
  // - choose one of these methods using `SetReconMethod(TString name)`
  // - if you specify none, a default method will be chosen
  reconMethodToTitle.insert({ "Ele",    "Electron method"        });
  reconMethodToTitle.insert({ "DA",     "Double Angle method"    });
  reconMethodToTitle.insert({ "JB",     "Jacquet-Blondel method" });
  reconMethodToTitle.insert({ "Mixed",  "Mixed method"           });
  reconMethodToTitle.insert({ "Sigma",  "Sigma method"           });
  reconMethodToTitle.insert({ "eSigma", "eSigma method"          });

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
// add a group of files
void Analysis::AddFileGroup(
    std::vector<std::string> fileNames,
    std::vector<Long64_t> entries,
    Double_t xs,
    Double_t Q2min,
    Double_t Q2max
    )
{
  // print
  fmt::print("-> AddFileGroup:  crossSection={},  Q2range = {} to {}\n",
      xs,
      Q2min,
      Q2max>Q2min ? std::to_string(Q2max) : "inf"
      );
  for (auto fileName : fileNames) fmt::print("   - {}\n", fileName);

  // fill vectors
  infiles.push_back(fileNames);
  Q2xsecs.push_back(xs);
  Q2mins.push_back(Q2min);
  Q2maxs.push_back(Q2max);
  inEntries.push_back(entries);

  // get total entries
  Long64_t totalEntry = 0;
  for (Long64_t entry : entries) totalEntry += entry;
  Q2entries.push_back(totalEntry);

  // check if the cross section for each group decreases; this is preferred
  // to make sure that `GetEventQ2Idx` behaves correctly
  for (std::size_t idx = 0; idx+1 < Q2xsecs.size(); ++idx) {
    if (Q2xsecs[idx] < Q2xsecs[idx+1] ) {
      cerr << "WARNING: Cross-sections should strictly decrease with stricter Q2min cuts; re-arrange your config file" << endl;
      PrintStdVector(Q2xsecs,"   cross sections");
    }
  }

  // lookup table, for tree number -> vector indices
  inLookup.clear();
  std::size_t total = 0;
  for (std::size_t idx = 0; idx < infiles.size(); ++idx) {
    total += infiles[idx].size();
    inLookup.resize(total, idx);
  }
}


// prepare for the analysis
//------------------------------------
void Analysis::Prepare() {

  // parse config file --------------------
  std::ifstream fin(infileName);
  std::string line;
  bool debugParser = false;
  fmt::print("[+++] PARSING CONFIG FILE {}\n",infileName);

  // vars
  Double_t xsec  = 0.0;
  Double_t Q2min = 1.0;
  Double_t Q2max = 0.0;
  std::vector<std::string> fileNames;
  bool readingListOfFiles = false;

  /* lambda to add a list of files to this analysis; wraps `Analysis::AddFileGroup`,
   * and checks the files beforehand
   */
  auto AddFiles = [this, &fileNames, &xsec, &Q2min, &Q2max] () {
    // get the number of entries, and check the file
    std::vector<Long64_t> entries;
    for(auto fileName : fileNames) {
      auto file = TFile::Open(fileName.c_str());
      if (file->IsZombie()) {
        fmt::print(stderr,"ERROR: Couldn't open input file '{}'\n",fileName);
        return;
      }
      TTree *tree = file->Get<TTree>("Delphes");                  // fastsim
      if (tree == nullptr) tree = file->Get<TTree>("events");     // ATHENA, EPIC
      if (tree == nullptr) tree = file->Get<TTree>("event_tree"); // ECCE
      if (tree == nullptr) {
        fmt::print(stderr,"ERROR: Couldn't find tree in file '{}'\n",fileName);
        entries.push_back(0);
      }
      else entries.push_back(tree->GetEntries());
      file->Close();
    }
    // add the file group
    AddFileGroup(fileNames, entries, xsec, Q2min, Q2max);
  };

  // loop over lines
  while (std::getline(fin, line)) {

    // chomp comments
    auto lineChomped = line.substr(0,line.find("#"));
    if(debugParser) fmt::print("\nlineChomped = \"{}\"\n",lineChomped);

    // each line will be parsed as one of the following:
    enum parseAs_enum { kNone, kSetting, kRootFile };
    int parseAs = kNone;

    // tokenize the line
    std::string token, bufferKey, bufferVal;
    std::stringstream lineStream(lineChomped);
    while (std::getline(lineStream, token, ' ')) {
      // classify `parseAs` and fill buffers
      if (token.rfind(":",0)==0) { // parse as setting name
        parseAs   = kSetting;
        bufferKey = token;
      }
      else if (parseAs==kSetting) { // parse as setting value
        if (token!="" && bufferVal=="") bufferVal = token;
      }
      else { // parse as root file name
        parseAs   = kRootFile;
        bufferKey = token;
      }
      // if(debugParser) fmt::print("  token = \"{}\"\n",token);
    } // tokenize

    // parse buffers
    if (parseAs==kSetting) {
      // if we were reading a list of files, we're done with this list; add the files and reset
      if(readingListOfFiles) {
        if(debugParser) fmt::print("-> new setting, add current list of `fileNames`, then reset\n");
        AddFiles();
        readingListOfFiles = false;
        xsec  = 0.0;
        Q2min = 1.0;
        Q2max = 0.0;
        fileNames.clear();
      }
      // parse setting value(s)
      if      (bufferKey==":eleBeamEn")         eleBeamEn         = std::stod(bufferVal);
      else if (bufferKey==":ionBeamEn")         ionBeamEn         = std::stod(bufferVal);
      else if (bufferKey==":crossingAngle")     crossingAngle     = std::stod(bufferVal);
      else if (bufferKey==":totalCrossSection") totalCrossSection = std::stod(bufferVal);
      else if (bufferKey==":q2min")             Q2min             = std::stod(bufferVal);
      else if (bufferKey==":q2max")             Q2max             = std::stod(bufferVal);
      else if (bufferKey==":crossSection")      xsec              = std::stod(bufferVal);
      else
        fmt::print(stderr,"WARNING: unkown setting \"{}\"\n",bufferKey);
      if(debugParser) fmt::print("  setting:  \"{}\" = \"{}\"\n",bufferKey,bufferVal);
    }
    else if (parseAs==kRootFile) {
      readingListOfFiles = true;
      fileNames.push_back(bufferKey);
      if(debugParser) fmt::print("  rootFile: \"{}\"\n",bufferKey);
    }

  } // loop over lines

  // add last list of files
  if(debugParser) fmt::print("-> EOF; add last list of `fileNames`\n");
  AddFiles();
  fmt::print("[+++] PARSING COMPLETE\n\n");

  if (infiles.empty()) {
    cerr << "ERROR: no input files have been specified" << endl;
    return;
  }

  // print configuration ----------------------------------
  fmt::print("{:=<50}\n","CONFIGURATION: ");
  fmt::print("{:>30} = {} GeV\n",  "eleBeamEn",         eleBeamEn);
  fmt::print("{:>30} = {} GeV\n",  "ionBeamEn",         ionBeamEn);
  fmt::print("{:>30} = {} mrad\n", "crossingAngle",     crossingAngle);
  fmt::print("{:>30} = {}\n",      "totalCrossSection", totalCrossSection);
  fmt::print("{:>30} = {}\n",      "reconMethod",       reconMethod);
  fmt::print("{:-<50}\n","");
  PrintStdVector(Q2mins,"Q2mins");
  PrintStdVector(Q2maxs,"Q2maxs");
  PrintStdVector(Q2xsecs,"Q2xsecs");
  PrintStdVector(Q2entries,"Q2entries");
  // PrintStdVector2D(inEntries,"inEntries");
  // PrintStdVector(inLookup,"inLookup");
  fmt::print("{:-<50}\n","");
  PrintStdMap(availableBinSchemes,"availableBinSchemes");
  fmt::print("{:=<50}\n","");

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

  // tell HD what values to associate for each BinScheme name
  // - these are the values that will be used for CutDef checks
  // - the lambda must return a Double_t
  /* DIS */
  HD->SetBinSchemeValue("x",     [this](){ return kin->x;                       });
  HD->SetBinSchemeValue("q2",    [this](){ return kin->Q2;                      });
  HD->SetBinSchemeValue("w",     [this](){ return kin->W;                       });
  HD->SetBinSchemeValue("y",     [this](){ return kin->y;                       });
  /* single hadron */
  HD->SetBinSchemeValue("p",     [this](){ return kin->pLab;                    });
  HD->SetBinSchemeValue("eta",   [this](){ return kin->etaLab;                  });
  HD->SetBinSchemeValue("pt",    [this](){ return kin->pT;                      });
  HD->SetBinSchemeValue("ptLab", [this](){ return kin->pTlab;                   });
  HD->SetBinSchemeValue("z",     [this](){ return kin->z;                       });
  HD->SetBinSchemeValue("qT",    [this](){ return kin->qT;                      });
  HD->SetBinSchemeValue("qTq",   [this](){ return kin->qT/TMath::Sqrt(kin->Q2); });
  HD->SetBinSchemeValue("mX",    [this](){ return kin->mX;                      });
  HD->SetBinSchemeValue("xF",    [this](){ return kin->xF;                      });
  HD->SetBinSchemeValue("phiH",  [this](){ return kin->phiH;                    });
  HD->SetBinSchemeValue("phiS",  [this](){ return kin->phiS;                    });
  HD->SetBinSchemeValue("tSpin", [this](){ return (Double_t)kin->tSpin;         });
  HD->SetBinSchemeValue("lSpin", [this](){ return (Double_t)kin->lSpin;         });
  /* jets */
#ifndef EXCLUDE_DELPHES
  HD->SetBinSchemeValue("ptJet", [this](){ return kin->pTjet;                   });
  HD->SetBinSchemeValue("zJet",  [this](){ return kin->zjet;                    });
#endif

  // some bin schemes values are checked here, instead of by CutDef checks ("External" cut type)
  // - the lambda must return a boolean
  HD->SetBinSchemeValueExternal("finalState", [this](Node *N){ return N->GetCut()->GetCutID() == finalStateID; });


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
        NBINS,1,3000,
        true,true
        );
    HS->DefineHist1D("Q2","Q2","GeV",NBINS,1.0,3000,true,true);
    HS->DefineHist1D("x","x","",NBINS,1e-3,1.0,true,true);
    HS->DefineHist1D("y","y","",NBINS,1e-3,1,true);
    HS->DefineHist1D("W","W","GeV",NBINS,0,50);
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
    HS->DefineHist1D("mX","m_{X}","GeV",NBINS,0,40);
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
    Double_t etabinsCoarse[] = {-4.0,-1.0,1.0,4.0};
    Double_t pbinsCoarse[] = {0.1,1,10,100};
    HS->DefineHist2D("etaVsPcoarse","p","#eta","GeV","",
        3, pbinsCoarse,
        3, etabinsCoarse,
        true,false
        );
    // -- depolarization
    HS->DefineHist2D("epsilonVsQ2", "Q^{2}", "#epsilon", "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 1.5, true, false);
    HS->DefineHist2D("depolAvsQ2",  "Q^{2}", "A",        "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 2.5, true, false);
    HS->DefineHist2D("depolBAvsQ2", "Q^{2}", "B/A",      "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 2.5, true, false);
    HS->DefineHist2D("depolCAvsQ2", "Q^{2}", "C/A",      "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 2.5, true, false);
    HS->DefineHist2D("depolVAvsQ2", "Q^{2}", "V/A",      "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 2.5, true, false);
    HS->DefineHist2D("depolWAvsQ2", "Q^{2}", "W/A",      "GeV^{2}", "", NBINS, 1, 3000, NBINS, 0, 2.5, true, false);
    // -- single-hadron cross sections
    //HS->DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
    HS->DefineHist1D("Q_xsec","Q","GeV",NBINS,1.0,3000,true,true); // log
    HS->Hist("Q_xsec")->SetMinimum(1e-10);
    // -- resolutions
    HS->DefineHist1D("x_Res","x-x_{true}","", NBINS, -0.5, 0.5);
    HS->DefineHist1D("y_Res","y-y_{true}","", NBINS, -0.2, 0.2);
    HS->DefineHist1D("Q2_Res","Q2-Q2_{true}","GeV^{2}", NBINS, -20, 20);
    HS->DefineHist1D("W_Res","W-W_{true}","GeV", NBINS, -20, 20);
    HS->DefineHist1D("Nu_Res","#nu-#nu_{true}","GeV", NBINS, -100, 100);
    HS->DefineHist1D("phiH_Res","#phi_{h}-#phi_{h}^{true}","", NBINS, -TMath::Pi(), TMath::Pi());
    HS->DefineHist1D("phiS_Res","#phi_{S}-#phi_{S}^{true}","", NBINS, -TMath::Pi(), TMath::Pi());
    HS->DefineHist1D("pT_Res","pT-pT^{true}","GeV", NBINS, -1.5, 1.5);
    HS->DefineHist1D("z_Res","z-z^{true}","", NBINS, -1.5, 1.5);
    HS->DefineHist1D("mX_Res","mX-mX^{true}","GeV", NBINS, -10, 10);
    HS->DefineHist1D("xF_Res","xF-xF^{true}","", NBINS, -1.5, 1.5);
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
    // -- jet kinematics
#ifndef EXCLUDE_DELPHES
    HS->DefineHist1D("pT_jet","jet p_{T}","GeV", NBINS, 1e-2, 50);
    HS->DefineHist1D("mT_jet","jet m_{T}","GeV", NBINS, 1e-2, 20);
    HS->DefineHist1D("z_jet","jet z","", NBINS,0, 1);
    HS->DefineHist1D("eta_jet","jet #eta_{lab}","", NBINS,-5,5);
    HS->DefineHist1D("qT_jet","jet q_{T}", "GeV", NBINS, 0, 10.0);
    HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
    HS->DefineHist1D("qTQ_jet","jet q_{T}/Q","", NBINS, 0, 3.0);
#endif
  });
  HD->ExecuteAndClearOps();


  // initialize total weights
  wTrackTotal = 0.;
  wJetTotal = 0.;
};

void Analysis::CalculateEventQ2Weights() {
  Q2weights.resize(Q2xsecs.size());
  entriesTot = 0;
  for (Long64_t entry : Q2entries) {
    entriesTot += entry;
  }
  cout << "Q2 weighting info:" << endl;
  for (Int_t idx = 0; idx < Q2xsecs.size(); ++idx) {
    // calculate total luminosity, and the luminosity that contains this Q2 range
    Double_t lumiTotal = Double_t(entriesTot)     / totalCrossSection;
    Double_t lumiThis  = Double_t(Q2entries[idx]) / Q2xsecs[idx];
    /*
    //// alternative `lumiThis`: try to correct for overlapping Q2 ranges; in practice this
    //// does not make much of a difference, in fact, the cross section looks slighly worse
    lumiThis = 0.;
    for (Int_t j = 0; j < Q2xsecs.size(); ++j) {
      // check if Q2 range `j` contains the Q2 range `idx`; if so, include its luminosity
      if (InQ2Range(Q2mins[idx], Q2mins[j], Q2maxs[j]) &&
          InQ2Range(Q2maxs[idx], Q2mins[j], Q2maxs[j], true))
      {
        lumiThis += Double_t(Q2entries[j]) / Q2xsecs[j];
      }
    }
    */
    // calculate the weight for this Q2 range
    Q2weights[idx] = lumiTotal / lumiThis;
    cout << "\tQ2 > "     << Q2mins[idx];
    if(Q2maxs[idx]>0) cout << " && Q2 < " << Q2maxs[idx];
    cout << ":" << endl;
    cout << "\t\tcount    = " << Q2entries[idx] << endl;
    cout << "\t\txsec     = " << Q2xsecs[idx] << endl;
    cout << "\t\tweight   = " << Q2weights[idx] << endl;
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

// get the `idx` for this value of `Q2`
Int_t Analysis::GetEventQ2Idx(Double_t Q2, Int_t guess) {
  Int_t idx = guess;
  // generally prefer the most restrictive Q2 range we have
  if (!InQ2Range(Q2,Q2mins[idx],Q2maxs[idx])) {
    // if Q2 is not in the expected range, try a less restrictive one
    do {
      idx -= 1;
    } while (idx >= 0 && !InQ2Range(Q2,Q2mins[idx],Q2maxs[idx]));
    // fmt::print("SCANNED DOWN to less restrictive Q2 range: Q2={}, guess={}, return idx={}\n",Q2,guess,idx);
    return idx;
  } else {
    // look for the most restrictive Q2 range that contains this Q2
    for (int k=idx+1; k<Q2mins.size(); k++) {
      if (InQ2Range(Q2,Q2mins[k],Q2maxs[k]))
        idx = k;
    }
    // if(idx!=guess) fmt::print("SCANNED UP to more restrictive Q2 range: Q2={}, guess={}, return idx={}\n",Q2,guess,idx);
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
  Double_t lumi = wTrackTotal/totalCrossSection; // [nb^-1]
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
  std::vector<Double_t> vec_wTrackTotal { wTrackTotal };
  std::vector<Double_t> vec_wJetTotal { wJetTotal };
  outFile->WriteObject(&Q2xsecs, "XsTotal");
  outFile->WriteObject(&vec_wTrackTotal, "WeightTotal");
  outFile->WriteObject(&vec_wJetTotal, "WeightJetTotal");

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
    binSchemes.insert({varname,B});
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


// FillHistos methods: check bins and fill associated histograms
// - checks which bins the track/jet/etc. falls in
// - fills the histograms in the associated Histos objects
//--------------------------------------------------------------------------
// tracks (single particles)
void Analysis::FillHistosTracks() {

  // check which bins to fill
  HD->CheckBins();

  // fill histograms, for activated bins only
  HD->Payload([this](Histos *H){
    // Full phase space.
    H->Hist4("full_xsec")->Fill(kin->x,kin->Q2,kin->pT,kin->z,wTrack);
    // DIS kinematics
    dynamic_cast<TH2*>(H->Hist("Q2vsX"))->Fill(kin->x,kin->Q2,wTrack);
    H->Hist("Q2")->Fill(kin->Q2,wTrack);
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
    dynamic_cast<TH2*>(H->Hist("etaVsPcoarse"))->Fill(kin->pLab,kin->etaLab,wTrack); 
    // depolarization
    dynamic_cast<TH2*>(H->Hist("epsilonVsQ2"))->Fill(kin->Q2,kin->epsilon,wTrack); 
    dynamic_cast<TH2*>(H->Hist("depolAvsQ2"))->Fill(kin->Q2,kin->depolA,wTrack); 
    dynamic_cast<TH2*>(H->Hist("depolBAvsQ2"))->Fill(kin->Q2,kin->depolP1,wTrack); 
    dynamic_cast<TH2*>(H->Hist("depolCAvsQ2"))->Fill(kin->Q2,kin->depolP2,wTrack); 
    dynamic_cast<TH2*>(H->Hist("depolVAvsQ2"))->Fill(kin->Q2,kin->depolP3,wTrack); 
    dynamic_cast<TH2*>(H->Hist("depolWAvsQ2"))->Fill(kin->Q2,kin->depolP4,wTrack); 
    // cross sections (divide by lumi after all events processed)
    H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),wTrack);
    // resolutions
    H->Hist("x_Res")->Fill( kin->x - kinTrue->x, wTrack );
    H->Hist("y_Res")->Fill( kin->y - kinTrue->y, wTrack );
    H->Hist("Q2_Res")->Fill( kin->Q2 - kinTrue->Q2, wTrack );
    H->Hist("W_Res")->Fill( kin->W - kinTrue->W, wTrack );
    H->Hist("Nu_Res")->Fill( kin->Nu - kinTrue->Nu, wTrack );
    H->Hist("phiH_Res")->Fill( Kinematics::AdjAngle(kin->phiH - kinTrue->phiH), wTrack );
    H->Hist("phiS_Res")->Fill( Kinematics::AdjAngle(kin->phiS - kinTrue->phiS), wTrack );
    H->Hist("pT_Res")->Fill( kin->pT - kinTrue->pT, wTrack );
    H->Hist("z_Res")->Fill( kin->z - kinTrue->z, wTrack );
    H->Hist("mX_Res")->Fill( kin->mX - kinTrue->mX, wTrack );
    H->Hist("xF_Res")->Fill( kin->xF - kinTrue->xF, wTrack );
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
#ifndef EXCLUDE_DELPHES
void Analysis::FillHistosJets() {

  // check which bins to fill
  HD->CheckBins();

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
#endif


// destructor
Analysis::~Analysis() {
};

