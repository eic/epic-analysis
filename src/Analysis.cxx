// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Connor Pecar, Duane Byer, Sanghwa Park, Matthew McEneaney, Brian Page

#include "Analysis.h"

ClassImp(Analysis)

using std::cout;
using std::cerr;
using std::endl;

// constructor
Analysis::Analysis(
  TString configFileName_,
  TString outfilePrefix_
)
  : configFileName(configFileName_)
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
  /* jets */
#ifndef EXCLUDE_DELPHES
  availableBinSchemes.insert({ "JetPT", "jet p_{T}" });
  availableBinSchemes.insert({ "JetZ",  "jet z"     });
  availableBinSchemes.insert({ "JetEta", "jet eta" });
  availableBinSchemes.insert({ "JetE", "jet energy" });
#endif

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

  // kinematics reconstruction methods
  // - choose one of these methods using `SetReconMethod(TString name)`
  // - if you specify none, a default method will be chosen
  reconMethodToTitle.insert({ "Ele",    "Electron method"        });
  reconMethodToTitle.insert({ "DA",     "Double Angle method"    });
  reconMethodToTitle.insert({ "JB",     "Jacquet-Blondel method" });
  reconMethodToTitle.insert({ "Mixed",  "Mixed method"           });
  reconMethodToTitle.insert({ "Sigma",  "Sigma method"           });
  reconMethodToTitle.insert({ "eSigma", "eSigma method"          });

  // output sets to include
  // - use these to turn on/off certain sets of variables
  // - the default settings are set here; override them at the macro level
  includeOutputSet.insert({ "inclusive",      true  }); // inclusive kinematics
  includeOutputSet.insert({ "1h",             true  }); // single hadron kinematics
  includeOutputSet.insert({ "jets",           false }); // jet kinematics
  includeOutputSet.insert({ "depolarization", false }); // depolarization factors & ratios

  // common settings defaults
  // - these settings can be set at the macro level
  verbose           = false;
  writeSidisTree   = false;
  writeHFSTree      = false;
  writeParticleTree = false;

  maxEvents    = 0;
  useBreitJets = false;
  errorCntMax  = 1000;
  jetAlg       = 0; // Default to kT Algorithm
  jetRad       = 0.8;
  jetMin       = 1.0; // Minimum Jet pT
  jetMatchDR   = 0.5; // Delta R between true and reco jet to be considered matched

  weightInclusive = std::make_unique<WeightsUniform>();
  weightTrack     = std::make_unique<WeightsUniform>();
  weightJet       = std::make_unique<WeightsUniform>();

  // miscellaneous
  infiles.clear();
  entriesTot = 0;
  errorCnt = 0;
  inputTreeName = "";
};


// input files
//------------------------------------
// add a group of files
void Analysis::AddFileGroup(
    std::vector<std::string> fileNames,
    Long64_t totalEntries,
    Double_t xs,
    Double_t Q2min,
    Double_t Q2max,
    Double_t manualWeight
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
  Q2entries.push_back(totalEntries);
  if (manualWeight!=0)
      Q2weights.push_back(manualWeight);
    
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
  std::ifstream fin(configFileName);
  std::string line;
  bool debugParser = false;
  fmt::print("[+++] PARSING CONFIG FILE {}\n",configFileName);

  // vars
  Double_t xsec;
  Double_t Q2min;
  Double_t Q2max;
  Double_t manualWeight;
  Long64_t numEvents;
  std::vector<std::string> fileNames;
  bool readingListOfFiles;
  bool endGroupCalled;

  // reset vars
  auto ResetVars = [&] () {
    xsec      = totalCrossSection;
    Q2min     = 1.0;
    Q2max     = 0.0;
    manualWeight    = 0.0;
    numEvents = -1;
    fileNames.clear();
    readingListOfFiles = false;
    endGroupCalled     = false;
  };
  ResetVars();

  /* lambda to add a list of files to this analysis; wraps `Analysis::AddFileGroup`,
   * and checks the files beforehand
   */
  auto AddFiles = [this, &fileNames, &xsec, &Q2min, &Q2max, &numEvents, &manualWeight] () {
    Long64_t entries = 0;
    if(numEvents<0) {
      // get the number of entries, and check the file
      for(auto fileName : fileNames) {
        auto file = TFile::Open(fileName.c_str());
        if (file==nullptr || file->IsZombie()) {
          fmt::print(stderr,"ERROR: Couldn't open input file '{}'\n",fileName);
          return;
        }
        TTree *tree = nullptr;
        if(inputTreeName == "") {
          // figure out which tree we are analyzing
          std::vector<std::string> inputTreeNameOpts = {
            "Delphes",   // fastsim
            "events",    // ePIC, ATHENA
            "event_tree" // ECCE
          };
          for (auto inputTreeNameOpt : inputTreeNameOpts) {
            tree = file->Get<TTree>(inputTreeNameOpt.c_str());
            if (tree != nullptr) {
              inputTreeName = inputTreeNameOpt;
              break;
            }
          }
          if (tree == nullptr)
            fmt::print(stderr,"ERROR: Couldn't find any known tree in file '{}'\n",fileName);
        }
        else {
          tree = file->Get<TTree>(inputTreeName.c_str());
          if (tree == nullptr)
            fmt::print(stderr,"ERROR: Couldn't find tree '{}' in file '{}'\n",inputTreeName,fileName);
        }
        if (tree != nullptr)
          entries += tree->GetEntries();
        file->Close();
        delete file;
      }
    }
    else entries = numEvents;
    // add the file group
    AddFileGroup(fileNames, entries, xsec, Q2min, Q2max, manualWeight);
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
      if(readingListOfFiles || bufferKey==":endGroup") {
        if(debugParser) fmt::print("-> new setting, add current list of `fileNames`, then reset\n");
        AddFiles();
        ResetVars();
      }
      // parse setting value(s)
      if(debugParser) fmt::print("  parse as setting: bufferKey='{}' bufferVal='{}'\n",bufferKey,bufferVal);
      if(bufferVal=="" && bufferKey!=":endGroup")
        fmt::print(stderr,"ERROR: setting '{}' has no associated value\n",bufferKey);
      else if (bufferKey==":eleBeamEn")         eleBeamEn         = std::stod(bufferVal);
      else if (bufferKey==":ionBeamEn")         ionBeamEn         = std::stod(bufferVal);
      else if (bufferKey==":crossingAngle")     crossingAngle     = std::stod(bufferVal);
      else if (bufferKey==":totalCrossSection") totalCrossSection = std::stod(bufferVal);
      else if (bufferKey==":q2min")             Q2min             = std::stod(bufferVal);
      else if (bufferKey==":q2max")             Q2max             = std::stod(bufferVal);
      else if (bufferKey==":crossSection")      xsec              = std::stod(bufferVal);
      else if (bufferKey==":Weight")            manualWeight      = std::stod(bufferVal);
      else if (bufferKey==":numEvents")         numEvents         = std::stoll(bufferVal);
      else if (bufferKey==":endGroup")          endGroupCalled    = true;
      else
        fmt::print(stderr,"ERROR: unknown setting \"{}\"\n",bufferKey);
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
  if(!endGroupCalled) AddFiles();
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
  // PrintStdVector(inLookup,"inLookup");
  // fmt::print("{:-<50}\n",""); PrintStdMap(availableBinSchemes,"availableBinSchemes");
  fmt::print("{:-<50}\n",""); PrintStdMap(includeOutputSet,"includeOutputSet");
  fmt::print("{:=<50}\n","");

  // set output file name
  outfileName = "out/"+outfilePrefix+".root";

  // open output file
  cout << "-- output file: " << outfileName << endl;
  outFile = new TFile(outfileName,"RECREATE");

  // instantiate shared objects
  kin        = std::make_shared<Kinematics>(eleBeamEn,ionBeamEn,crossingAngle);
  kinTrue    = std::make_shared<Kinematics>(eleBeamEn, ionBeamEn, crossingAngle);
#ifndef EXCLUDE_DELPHES
  kinJet     = std::make_shared<KinematicsJets>(eleBeamEn, ionBeamEn, crossingAngle);
  kinJetTrue = std::make_shared<KinematicsJets>(eleBeamEn, ionBeamEn, crossingAngle);
#endif
  ST         = std::make_unique<SidisTree>("tree",kin,kinTrue);
  HFST       = std::make_unique<HFSTree>("hfstree",kin,kinTrue);
  PT         = std::make_unique<ParticleTree>("ptree");

  // if including jets, define a `jet` final state
#ifndef EXCLUDE_DELPHES
  if(includeOutputSet["jets"]) {
    finalStateToTitle.insert({ "jet", "jets" });
    AddFinalState("jet");
  }
#endif

  // determine if only including `inclusive` output set
  includeOutputSet.insert({ "inclusive_only",
      includeOutputSet["inclusive"]
      && !includeOutputSet["1h"]
      && !includeOutputSet["jets"]
      });

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
  HD = std::make_shared<HistosDAG>();
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
  HD->SetBinSchemeValue("JetPT", [this](){ return kinJet->pTjet;                });
  HD->SetBinSchemeValue("JetZ",  [this](){ return kinJet->zjet;                 });
  HD->SetBinSchemeValue("JetEta", [this](){ return kinJet->etajet;              });
  HD->SetBinSchemeValue("JetE", [this](){ return kinJet->ejet;                  });
#endif

  // some bin schemes values are checked here, instead of by CutDef checks ("External" cut type)
  // - the lambda must return a boolean
  HD->SetBinSchemeValueExternal("finalState", [this](Node *N){ return N->GetCut()->GetCutID() == finalStateID; });


  // DEFINE HISTOGRAMS ------------------------------------
  // - whether they are defined is controlled by `includeOutputSet` settings
  // - `Histos::Hist` calls should check for existence
  HD->Payload([this](Histos *HS){
    // -- inclusive kinematics
    if(includeOutputSet["inclusive"]) {
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
    }

    // -- single hadron kinematics
    if(includeOutputSet["1h"]) {
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
          true,false,true
          );
      Double_t etabinsCoarse[] = {-4.0,-1.0,1.0,4.0};
      Double_t pbinsCoarse[] = {0.1,1,10,100};
      HS->DefineHist2D("etaVsPcoarse","p","#eta","GeV","",
          3, pbinsCoarse,
          3, etabinsCoarse,
          true,false
          );
      // -- single-hadron cross sections
      //HS->DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
      HS->DefineHist1D("Q_xsec","Q","GeV",NBINS,1.0,3000,true,true); // log
      HS->Hist("Q_xsec")->SetMinimum(1e-10);
      // -- resolutions
      HS->DefineHist1D("x_Res",  "x-x_{true}/x_{true}",             "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("y_Res",  "y-y_{true}/y_{true}",             "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("Q2_Res", "Q^{2}-Q^{2}_{true}/Q^{2}_{true}", "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("W_Res",  "W-W_{true}/W_{true}",             "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("Nu_Res", "#nu-#nu_{true}/#nu_{true}",       "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("pT_Res", "p_{T}-p_{T}^{true}/p_{T}^{true}", "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("z_Res",  "z-z_{true}/z_{true}",             "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("mX_Res", "m_{X}-m_{X}^{true}/m_{X}^{true}", "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("xF_Res", "x_{F}-x_{F}^{true}/x_{F}^{true}", "", NBINS, -1.0, 1.0);
      HS->DefineHist1D("phiH_Res", "#phi_{h}-#phi_{h}^{true}", "", NBINS, -TMath::Pi(), TMath::Pi()); // absolute resolution
      HS->DefineHist1D("phiS_Res", "#phi_{S}-#phi_{S}^{true}", "", NBINS, -TMath::Pi(), TMath::Pi()); // absolute resolution
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
    }

    // -- jet kinematics
#ifndef EXCLUDE_DELPHES
    if(includeOutputSet["jets"]) {
      HS->DefineHist1D("JetPT",    "jet p_{T}",      "GeV", NBINS, 1e-2, 50,  true, true);
      HS->DefineHist1D("JetMT",    "jet m_{T}",      "GeV", NBINS, 1e-2, 20,  true, true);
      HS->DefineHist1D("JetZ",     "jet z",          "",    NBINS, 0,    1);
      HS->DefineHist1D("JetEta",   "jet #eta_{lab}", "",    NBINS, -5,   5);
      HS->DefineHist1D("JetQT",    "jet q_{T}",      "GeV", NBINS, 0,    10.0);
      HS->DefineHist1D("JetJperp", "j_{#perp}",      "GeV", NBINS, 1e-2, 3.0, true, true);
      HS->DefineHist1D("JetQTQ",   "jet q_{T}/Q",    "",    NBINS, 0,    3.0);
      HS->DefineHist1D("JetE","jet Energy","GeV", NBINS, 1e-2, 100);
      HS->DefineHist1D("JetM","jet Mass","GeV", NBINS, 1e-2, 20);

      HS->DefineHist2D("JetPTVsEta","Eta","pT","GeV","",100,-5,5,100,0,50,false,false);

      // Resolution Plos
      HS->DefineHist1D("JetDeltaR","deltaR between matched reco and truth jet","",1000,-2.,10.);
 
      HS->DefineHist2D("JetPTTrueVsReco","Reco pT","True pT","GeV","GeV",100,0.,50.,100,0.,50.,false,false);
      HS->DefineHist2D("JetETrueVsReco","Reco E","True E","GeV","GeV",100,0.,100.,100,0.,100.,false,false);
 
      HS->DefineHist2D("JetResPTVsTruePT","True pT","(Reco-True)/True pT","GeV","",100,0.,50.,10000,-10.,10.,false,false);
      HS->DefineHist2D("JetResEVsTrueE","True E","(Reco-True)/True E","GeV","",100,0.,100.,10000,-10.,10.,false,false);
      
      HS->DefineHist2D("JetResPTVsRecoPT","Reco pT","(Reco-True)/True pT","GeV","",100,0.,50.,10000,-10.,10.,false,false);
      HS->DefineHist2D("JetResEVsRecoE","Reco E","(Reco-True)/True E","GeV","",100,0.,100.,10000,-10.,10.,false,false);
    }
#endif

    // -- depolarization
    if(includeOutputSet["depolarization"]) {
      std::map<TString,TString> depols = {
        { "epsilon", "#epsilon" },
        { "depolA",  "A"        },
        { "depolB",  "B"        },
        { "depolV",  "V"        },
        { "depolBA", "B/A"      },
        { "depolCA", "C/A"      },
        { "depolVA", "V/A"      },
        { "depolWA", "W/A"      }
      };
      for(auto [name,title] : depols) {
        HS->DefineHist2D(name+"vsQ2", "Q^{2}", title, "GeV^{2}", "", NBINS, 1,    3000, NBINS, 0, 2.5, true, false);
        HS->DefineHist2D(name+"vsY",  "y",     title, "",        "", NBINS, 5e-3, 1,    NBINS, 0, 2.5, true, false);
        HS->DefineHist3D(name+"vsQ2vsX",
            "x",   "Q^{2}",   title,
            "",    "GeV^{2}", "",
            NBINS, 1e-3,      1,
            NBINS, 1,         3000,
            NBINS, 0,         2.5,
            true,  true,      false
            );
      }
    }

  });
  HD->ExecuteAndClearOps();


  // initialize total weights
  wInclusiveTotal = 0.;
  wTrackTotal     = 0.;
  wJetTotal       = 0.;
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
    //// alternative `lumiThis`: try to correct for overlapping Q2 ranges
    lumiThis = 0.; // reset
    for (Int_t j = 0; j < Q2xsecs.size(); ++j) {
      // check if Q2 range `j` contains the Q2 range `idx`; if so, include its luminosity
      if (InQ2Range(Q2mins[idx], Q2mins[j], Q2maxs[j]) &&
          InQ2Range(Q2maxs[idx], Q2mins[j], Q2maxs[j], true))
      {
        lumiThis += Double_t(Q2entries[j]) / Q2xsecs[j];
      }
    }
    // calculate the weight for this Q2 range
    if(Q2weights[idx]==0)
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

// get the `idx` for this value of `Q2`; tries to find the most restrictive Q2 range
// that contains the given `Q2` value, in order to assign the correct weight
Int_t Analysis::GetEventQ2Idx(Double_t Q2, Int_t guess) {
  // return guess; // just use the Q2 range, according to which ROOT file the event came from
  Int_t idx = -1;
  for (int k=0; k<Q2mins.size(); k++) {
    if (InQ2Range(Q2,Q2mins[k],Q2maxs[k])) idx = k;
  }
  if (idx<0) {
    // fmt::print(stderr,"WARNING: Q2={} not in any Q2 range\n",Q2); // usually just below the smallest Q2 min
    idx = 0; // assume the least-restrictive range
  }
  else if(idx<guess){
      // When crossingAngle!=0, some event can be reconstructed with trueQ2 < Q2min of the Monte Carlo.
      // If so, set the idx=guess, i.e. just assume the event was generated in its file's Q2range
      idx = guess; 
  }
  return idx;
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
    auto h_Q2vsX       = H->Hist("Q2vsX",true);
    auto h_Q_xsec      = H->Hist("Q_xsec",true);
    auto h_Q2vsXpurity = H->Hist("Q2vsXpurity",true);
    if(h_Q2vsX) {
      cout << H->GetSetTitle() << " ::: "
           << h_Q2vsX->GetEntries()
           << endl;
    }
    // calculate cross sections
    if(h_Q_xsec) h_Q_xsec->Scale(1./lumi); // TODO: generalize (`if (name contains "xsec") ...`)
    // divide resolution plots by true counts per x-Q2 bin
    if(h_Q2vsXpurity) {
      H->Hist("Q2vsXpurity")->Divide(H->Hist("Q2vsXtrue"));
      H->Hist("Q2vsX_zres")->Divide(H->Hist("Q2vsXtrue"));
      H->Hist("Q2vsX_pTres")->Divide(H->Hist("Q2vsXtrue"));
      H->Hist("Q2vsX_phiHres")->Divide(H->Hist("Q2vsXtrue"));        
    }
  });
  HD->ExecuteAndClearOps();

  // write histograms
  cout << sep << endl;
  cout << "writing ROOT file..." << endl;
  outFile->cd();
  if(writeSidisTree)   ST->WriteTree();
  if(writeHFSTree)      HFST->WriteTree();
  if(writeParticleTree) PT->WriteTree();
  HD->Payload([this](Histos *H){ H->WriteHists(outFile); }); HD->ExecuteAndClearOps();
  HD->Payload([this](Histos *H){ H->Write(); }); HD->ExecuteAndClearOps();
  std::vector<Double_t> vec_wInclusiveTotal { wInclusiveTotal };
  std::vector<Double_t> vec_wTrackTotal     { wTrackTotal     };
  std::vector<Double_t> vec_wJetTotal       { wJetTotal       };
  outFile->WriteObject(&Q2xsecs, "XsTotal");
  outFile->WriteObject(&vec_wInclusiveTotal, "WeightInclusiveTotal");
  outFile->WriteObject(&vec_wTrackTotal,     "WeightTrackTotal");
  outFile->WriteObject(&vec_wJetTotal,       "WeightJetTotal");
  outFile->WriteObject(&total_events,    "TotalEvents");
    
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
  fmt::print("AddFinalState: name='{}'\n               title='{}'\n",finalStateN,finalStateT);
};


// FillHistos methods: check bins and fill associated histograms
// - checks which bins the track/jet/etc. falls in
// - fills the histograms in the associated Histos objects
//--------------------------------------------------------------------------

// fill histograms (called by other specific FillHistos* methods)
void Analysis::FillHistos(std::function<void(Histos*)> fill_payload) {
  // check which bins to fill, and activate the ones for which all defined cuts pass
  // (activates their corresponding `HistosDAG` nodes)
  HD->CheckBins();
  // fill histograms, for activated bins only
  HD->Payload(fill_payload);
  // execute the payload
  // - save time and don't call `ClearOps` (next loop will overwrite lambda)
  // - called with `activeNodesOnly==true` since we only want to fill bins associated
  //   with this track
  HD->ExecuteOps(true);
};

// fill inclusive histograms
void Analysis::FillHistosInclusive(Double_t wgt) {
  auto fill_payload = [this,wgt] (Histos *H) {
    H->FillHist2D("Q2vsX", kin->x,  kin->Q2, wgt);
    H->FillHist1D("Q2",    kin->Q2, wgt);
    H->FillHist1D("x",     kin->x,  wgt);
    H->FillHist1D("W",     kin->W,  wgt);
    H->FillHist1D("y",     kin->y,  wgt);
  };
  FillHistos(fill_payload);
};

// fill 1h (single-hadron) histograms
void Analysis::FillHistos1h(Double_t wgt) {
  auto fill_payload = [this,wgt] (Histos *H) {
    // Full phase space.
    H->FillHist4D("full_xsec", kin->x, kin->Q2, kin->pT, kin->z, wgt);
    // hadron 4-momentum
    H->FillHist1D("pLab",   kin->pLab,   wgt);
    H->FillHist1D("pTlab",  kin->pTlab,  wgt);
    H->FillHist1D("etaLab", kin->etaLab, wgt);
    H->FillHist1D("phiLab", kin->phiLab, wgt);
    // hadron kinematics
    H->FillHist1D("z",  kin->z,  wgt);
    H->FillHist1D("pT", kin->pT, wgt);
    H->FillHist1D("qT", kin->qT, wgt);
    if(kin->Q2!=0) H->FillHist1D("qTq",kin->qT/TMath::Sqrt(kin->Q2),wgt);
    H->FillHist1D("mX",         kin->mX,   wgt);
    H->FillHist1D("phiH",       kin->phiH, wgt);
    H->FillHist1D("phiS",       kin->phiS, wgt);
    H->FillHist2D("phiHvsPhiS", kin->phiS, kin->phiH, wgt);
    H->FillHist1D("phiSivers",  Kinematics::AdjAngle(kin->phiH - kin->phiS), wgt);
    H->FillHist1D("phiCollins", Kinematics::AdjAngle(kin->phiH + kin->phiS), wgt);
    H->FillHist2D("etaVsP",       kin->pLab, kin->etaLab, wgt); // TODO: lab-frame p, or some other frame?
    H->FillHist2D("etaVsPcoarse", kin->pLab, kin->etaLab, wgt);
    // depolarization
    if(includeOutputSet["depolarization"]) { // not necessary, but done for performance
      std::map<TString,Double_t> depols = {
        { "epsilon", kin->epsilon },
        { "depolA",  kin->depolA  },
        { "depolB",  kin->depolB  },
        { "depolV",  kin->depolV  },
        { "depolBA", kin->depolP1 },
        { "depolCA", kin->depolP2 },
        { "depolVA", kin->depolP3 },
        { "depolWA", kin->depolP4 }
      };
      for(auto [name,val] : depols) {
        H->FillHist2D(name+"vsQ2",    kin->Q2, val,     wgt);
        H->FillHist2D(name+"vsY",     kin->y,  val,     wgt);
        H->FillHist3D(name+"vsQ2vsX", kin->x,  kin->Q2, val, wgt);
      }
    }
    // cross sections (divide by lumi after all events processed)
    H->FillHist1D("Q_xsec", TMath::Sqrt(kin->Q2), wgt);
    // resolutions
    auto fillRelativeResolution = [&H,&wgt] (auto name, auto rec, auto gen) { if(gen!=0) H->FillHist1D(name, (rec-gen)/gen, wgt); };
    fillRelativeResolution("x_Res",  kin->x,  kinTrue->x);
    fillRelativeResolution("y_Res",  kin->y,  kinTrue->y);
    fillRelativeResolution("Q2_Res", kin->Q2, kinTrue->Q2);
    fillRelativeResolution("W_Res",  kin->W,  kinTrue->W);
    fillRelativeResolution("Nu_Res", kin->Nu, kinTrue->Nu);
    fillRelativeResolution("pT_Res", kin->pT, kinTrue->pT);
    fillRelativeResolution("z_Res",  kin->z,  kinTrue->z);
    fillRelativeResolution("mX_Res", kin->mX, kinTrue->mX);
    fillRelativeResolution("xF_Res", kin->xF, kinTrue->xF);
    H->FillHist1D("phiH_Res",  Kinematics::AdjAngle(kin->phiH - kinTrue->phiH), wgt );
    H->FillHist1D("phiS_Res",  Kinematics::AdjAngle(kin->phiS - kinTrue->phiS), wgt );
    H->FillHist2D("Q2vsXtrue", kinTrue->x,            kinTrue->Q2, wgt);
    if(kinTrue->z!=0)
      H->FillHist2D("Q2vsX_zres", kinTrue->x, kinTrue->Q2, wgt*( fabs(kinTrue->z - kin->z)/(kinTrue->z) ) );
    if(kinTrue->pT!=0)
      H->FillHist2D("Q2vsX_pTres", kinTrue->x, kinTrue->Q2, wgt*( fabs(kinTrue->pT - kin->pT)/(kinTrue->pT) ) );
    H->FillHist2D("Q2vsX_phiHres", kinTrue->x, kinTrue->Q2, wgt*( fabs(Kinematics::AdjAngle(kinTrue->phiH - kin->phiH) ) ) );

    auto htrue = H->Hist("Q2vsXtrue",true);
    if(htrue!=nullptr) {
      if( htrue->FindBin(kinTrue->x,kinTrue->Q2) == htrue->FindBin(kin->x,kin->Q2) )
        H->FillHist2D("Q2vsXpurity", kin->x, kin->Q2, wgt);
    }
    // -- reconstructed vs. generated
    H->FillHist2D("x_RvG",    kinTrue->x,    kin->x,    wgt);
    H->FillHist2D("phiH_RvG", kinTrue->phiH, kin->phiH, wgt);
    H->FillHist2D("phiS_RvG", kinTrue->phiS, kin->phiS, wgt);
  };
  FillHistos(fill_payload);
};

// fill jet histograms
void Analysis::FillHistosJets(Double_t wgt) {
#ifndef EXCLUDE_DELPHES
  auto fill_payload = [this,wgt] (Histos *H) {
    // jet kinematics
    H->FillHist1D("JetPT",  kinJet->pTjet,  wgt);
    H->FillHist1D("JetMT",  kinJet->mTjet,  wgt);
    H->FillHist1D("JetZ",   kinJet->zjet,   wgt);
    H->FillHist1D("JetEta", kinJet->etajet, wgt);
    H->FillHist1D("JetQT",  kinJet->qTjet,  wgt);
    H->FillHist1D("JetE",   kinJet->ejet,   wgt);
    H->FillHist1D("JetM",   kinJet->mjet,   wgt);

    H->FillHist2D("JetPTVsEta", kinJet->etajet, kinJet->pTjet, wgt);

    // Fill Resolution Histos
    H->FillHist1D("JetDeltaR", kinJet->deltaRjet, wgt);

    H->FillHist2D("JetPTTrueVsReco", kinJet->pTjet, kinJet->pTmtjet, wgt);
    H->FillHist2D("JetETrueVsReco", kinJet->ejet, kinJet->emtjet, wgt);

    H->FillHist2D("JetResPTVsTruePT", kinJet->pTmtjet, (kinJet->pTjet - kinJet->pTmtjet)/(kinJet->pTmtjet), wgt);
    H->FillHist2D("JetResEVsTrueE", kinJet->emtjet, (kinJet->ejet - kinJet->emtjet)/(kinJet->emtjet), wgt);

    H->FillHist2D("JetResPTVsRecoPT", kinJet->pTjet, (kinJet->pTjet - kinJet->pTmtjet)/(kinJet->pTmtjet), wgt);
    H->FillHist2D("JetResEVsRecoE", kinJet->ejet, (kinJet->ejet - kinJet->emtjet)/(kinJet->emtjet), wgt);

    if(kinJet->Q2!=0) H->FillHist1D("JetQTQ", kinJet->qTjet/sqrt(kinJet->Q2), wgt);
    for(int j = 0; j < kinJet->jperp.size(); j++) {
      H->FillHist1D("JetJperp", kinJet->jperp[j], wgt);
    };
  };
  FillHistos(fill_payload);
#endif
};


// print an error; if more than `errorCntMax` errors are printed, printing is suppressed
void Analysis::ErrorPrint(std::string message) {
  errorCnt++;
  if(errorCnt <= errorCntMax) fmt::print(stderr,"{}\n",message);
  if(errorCnt == errorCntMax) fmt::print(stderr,"... {} errors printed; suppressing the rest ...\n",errorCnt);
}

// destructor
Analysis::~Analysis() {
  for(auto it : binSchemes)
    if(it.second) delete it.second;
  binSchemes.clear();
}
