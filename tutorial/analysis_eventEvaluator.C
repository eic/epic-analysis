// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

/* full simulation (EventEvaluator/Fun4all) usage
 * - note the similarity of the macro to the fast simulation
 * - you only need to swap `AnalysisDelphes` with `AnalysisEcce` to switch
 *   between fast and full simulations
 * - some settings are specific to the full simulations, e.g. electron
 *   energy threshold
 * - this tutorial accesses files on S3:
 *   - alternatively, use your preferred method to run (download, mirror, etc.)
 *   - for S3, you must know the username and password, and have them in your environment:
 *     - `export S3_ACCESS_KEY=<login>`
 *     - `export S3_SECRET_KEY=<password>`
 */
void analysis_eventEvaluator(
    TString configFile="tutorial/ecce.config",      // input config file
    TString outfilePrefix="tutorial.eventEvaluator" // output filename prefix
) {

  // setup analysis ========================================
  // - define `AnalysisEcce` instead of `AnalysisDelphes`
  AnalysisEcce *A = new AnalysisEcce(configFile, outfilePrefix);

  A->maxEvents = 300000; // use this to limit the number of events
  A->writeSidisTree = true;

  // set reconstruction method and final states =============================
  // - see `Analysis` constructor for methods (or other tutorials)
  A->SetReconMethod("Ele");

  // decide which output sets to include ===================
  // - by default, only single-hadron data are included in the output
  // - see `src/Analysis.cxx` for all available settings
  //A->includeOutputSet["jets"] = true;
  //A->includeOutputSet["1h"] = false;
  //A->includeOutputSet["inclusive"] = false;
  //A->includeOutputSet["depolarization"] = true;

  A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");


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
  // - see other tutorials for guidance
  // - see `Analysis` constructor for available bin variables
  A->AddBinScheme("q2");
  A->AddBinScheme("x");

  // 3x2 grid of (x,Q2) bins, equal width in logarithmic scale
  A->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );



  // perform the analysis ==================================
  A->Execute();

  // for reference, here is a print out of HistosDAG
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");

}
