// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

/* run in a grid of (x,Q2,eta,p) 4D bins
 */
void analysis_coverage(
    TString configFile="tutorial/delphes.config", // input config file
    TString outfilePrefix="coverage"              // output filename prefix
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(configFile, outfilePrefix);

  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state


  // define cuts ====================================
  // - cuts are defined the same way as bins are defined; be mindful
  //   of what bins you are defining vs. what cuts you are defining.
  //   For example, if you define a Q2 minimum cut, and you also define
  //   Q2 bins below, you may be creating more bins than you actually
  //   need, since the Q2 minimum cut actually defines another bin;
  //   in this case, your Q2 bins effectively define a Q2 minimum.
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================
  /* - we define 4-D bins in (Q2,x,eta,p)
   * - `BuildBins` is used for each variable, 2 bins each
   * - a "Full" bin is added for each variable, where "Full"
   *   means no bin cut is applied; this will allow us to 
   *   "integrate" over variables we don't care about
   * - for example, make "eta vs. p" plots in bins of
   *   (x,Q2): since eta and p bins are also defined, we would
   *   only want to make this plot when in the "Full" eta
   *   and "Full" p bins, with a subloop through all (x,Q2) bins
   * - to achieve this behavior, we will use "ConditionalControl"
   *   in the postprocessor macro
   * - there are 3x3x3x3=81 bins
   */
  A->AddBinScheme("q2");
  A->BinScheme("q2")->BuildBin("Full");
  A->BinScheme("q2")->BuildBins( 2, 1, 100, true );

  A->AddBinScheme("x");
  A->BinScheme("x")->BuildBin("Full");
  A->BinScheme("x")->BuildBins( 2, 0.01, 1, true );

  A->AddBinScheme("eta");
  A->BinScheme("eta")->BuildBin("Full");
  A->BinScheme("eta")->BuildBins( 2, -0.5, 4.0 );

  A->AddBinScheme("p");
  A->BinScheme("p")->BuildBin("Full");
  A->BinScheme("p")->BuildBins( 2, 0.01, 10, true );

  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
