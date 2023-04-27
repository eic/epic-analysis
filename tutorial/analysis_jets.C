// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Brian Page

R__LOAD_LIBRARY(EpicAnalysis)

/* template analysis macro
 * - runs an analysis and nothing more (e.g., it does not define bins)
 * - useful if you want to use your own analysis code
 */
void analysis_jets(
    TString configFile="tutorial/delphes.config", // input config file
    TString outfilePrefix="tutorial.template"     // output filename prefix
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(configFile, outfilePrefix);
  // use a different `Analysis`-derived class to handle different data streams:
  // AnalysisEpic *A = new AnalysisEpic(configFile, outfilePrefix);     // ePIC Single Software Stack
  // AnalysisEcce *A = new AnalysisEcce(configFile, outfilePrefix);     // ECCE Fun4all + EventEvaluator
  // AnalysisAthena *A = new AnalysisAthena(configFile, outfilePrefix); // ATHENA DD4hep + Juggler

  //A->maxEvents = 10000; // use this to limit the number of events
  A->writeSidisTree = false; // (don't) write SidisTree (for one bin)

  // set reconstruction method =============================
  // - see `Analysis` constructor for methods
  A->SetReconMethod("Ele");

  // Include Jet output set and disable others ===================
  A->includeOutputSet["jets"] = true;
  // - additional example settings; see `src/Analysis.cxx` for more
  A->includeOutputSet["1h"] = false;
  A->includeOutputSet["inclusive"] = false;
  A->includeOutputSet["depolarization"] = false;

  // Define Jet Parameters ==================================
  // jetAlg controls the algorithm used (jetAlg=0 -> kT (default), jetAlg=1 -> Cambridge, jetAlg=2 -> anti_kT
  // Any other value will result in default being used
  A->jetAlg = 2; // Set Jet Algorithm
  A->jetRad = 1.0; // Set Jet Recombination Parameter
  A->jetMin = 1.0; // Set Minimum pT for Found Jet
  A->jetMatchDR = 0.5; // Set Minimum Delta R for Reco - True Jet Matching

  // Define Event Level Cuts ====================================
  // Keep loose for jets, only define a 'y' cut now
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  
  // Define Jet Level Cuts
  A->AddBinScheme("JetPT"); A->BinScheme("JetPT")->BuildBin("Min",5.0);
  A->AddBinScheme("JetEta");
  A->BinScheme("JetEta")->BuildBin("Range",-5.0,5.0);
  A->BinScheme("JetEta")->BuildBin("Range",-5.0,-1.0);
  A->BinScheme("JetEta")->BuildBin("Range",-1.0,1.0);
  A->BinScheme("JetEta")->BuildBin("Range",1.0,5.0);

  // set binning scheme ====================================
  // - see `Analysis` constructor for available bin variables
  /* do nothing -> single bin histograms */

  // perform the analysis ==================================
  A->Execute();
};
