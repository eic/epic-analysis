// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

/* template analysis macro
 * - runs an analysis and nothing more (e.g., it does not define bins)
 * - useful if you want to use your own analysis code
 */
void analysis_template(
    TString configFile="tutorial/delphes.config", // input config file
    TString outfilePrefix="tutorial.template"     // output filename prefix
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(configFile, outfilePrefix);
  // use a different `Analysis`-derived class to handle different data streams:
  //AnalysisEpic   *A = new AnalysisEpic(configFile,   outfilePrefix); // ePIC Single Software Stack
  //AnalysisEcce   *A = new AnalysisEcce(configFile,   outfilePrefix); // ECCE Fun4all + EventEvaluator
  //AnalysisAthena *A = new AnalysisAthena(configFile, outfilePrefix); // ATHENA DD4hep + Juggler

  //A->maxEvents = 10000; // use this to limit the number of events
  A->writeSidisTree = true; // write SidisTree (for one bin)

  // set reconstruction method =============================
  /* - currently we have:
   *   - "Ele"    => Electron method
   *   - "DA"     => Double Angle method
   *   - "JB"     => Jacquet-Blondel method
   *   - "Mixed"  => Mixed method
   *   - "Sigma"  => Sigma method
   *   - "eSigma" => eSigma method
   * - see `Analysis` constructor in `src/Analysis.cxx` for the most up-to-date list
   */
  A->SetReconMethod("Ele");

  // decide which output sets to include ===================
  // - by default, only single-hadron data are included
  // - to include jets, we just need to make sure they are included in the output
  A->includeOutputSet["jets"] = true;
  // - additional example settings; see `src/Analysis.cxx` for more
  //A->includeOutputSet["1h"] = false;            // exclude single-hadron data
  //A->includeOutputSet["inclusive"] = false;     // exclude inclusive data
  //A->includeOutputSet["depolarization"] = true; // include plots of depolarization factors

  // define cuts ====================================
  // - cuts are defined the same way as bins are defined; be mindful
  //   of what bins you are defining vs. what cuts you are defining.
  //   For example, if you define a Q2 minimum cut, and you also define
  //   Q2 bins below, you may be creating more bins than you actually
  //   need, since the Q2 minimum cut actually defines another bin;
  //   in this case, your Q2 bins effectively define a Q2 minimum.
  // - the cuts listed here are for single-hadrons only; similar cuts can be
  //   defined for jets
  A->AddBinScheme("w");     A->BinScheme("w")->BuildBin("Min",3.0);         // W > 3 GeV
  A->AddBinScheme("y");     A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");     A->BinScheme("z")->BuildBin("Range",0.2,0.9);   // 0.2 < z < 0.9
  A->AddBinScheme("xF");    A->BinScheme("xF")->BuildBin("Min",0.0);        // xF > 0
  A->AddBinScheme("ptLab"); A->BinScheme("ptLab")->BuildBin("Min",0.1);     // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  // - see `Analysis` constructor for available bin variables
  // - other tutorial macros define binning schemes
  /* do nothing -> single bin histograms */

  // final states =========================================
  // - define single-hadron final states; if you define none, default sets will be defined
  // - note: if you included jets above, there will also be a final state for jets, automatically defined
  // - see `src/Analysis.cxx` for more final states
  A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");

  // perform the analysis ==================================
  A->Execute();
};
