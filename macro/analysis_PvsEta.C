// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Connor Pecar

R__LOAD_LIBRARY(EpicAnalysis)

/* run in a grid of (x,Q2) 2D bins
 * - various ways to make a grid are demonstrated
 * - observe how the resulting histograms differ in each (x,Q2) bin
 */
void analysis_PvsEta(
    TString configFile="datarec/in.config", /* delphes tree(s) */
    TString outfilePrefix="coverage_pVsEtabins" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      configFile,
      outfilePrefix
      );

  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  A->AddFinalState("KpTrack"); // kaon final state


  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================

  /* TODO
   * - finer binning
   * - different sqrt(s) values
   * - eta vs. p in bins of (x,Q2)
   * - Q2 vs. x in bins of (eta,p)
   */

  A->AddBinScheme("p");
  A->AddBinScheme("eta");
  A->AddBinScheme("z");

  A->BinScheme("p")->BuildBin( "Range", 0, 2);
  A->BinScheme("p")->BuildBin( "Range", 2, 4);
  A->BinScheme("p")->BuildBin( "Range", 4, 10);
  A->BinScheme("p")->BuildBin( "Min", 10);

  A->BinScheme("eta")->BuildBin( "Range", 2.0, 3.5);
  A->BinScheme("eta")->BuildBin( "Range", 1.0, 2.0);
  A->BinScheme("eta")->BuildBin( "Range", -1.0, 1.0);
  A->BinScheme("eta")->BuildBin( "Range", -1.0, -2.0);
  A->BinScheme("eta")->BuildBin( "Range", -2.0, -3.5);

  A->BinScheme("z")->BuildBin("Range", 0.2, 0.4 );
  A->BinScheme("z")->BuildBin("Range", 0.4, 0.8 );

  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
