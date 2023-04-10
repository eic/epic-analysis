// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// run in Q2 bins, for two pT ranges
void analysis_qbins(
    TString configFile="tutorial/delphes.config", // input config file
    TString outfilePrefix="tutorial.qbins"        // output filename prefix
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
  A->AddBinScheme("q2");
  A->BinScheme("q2")->BuildBins( 10, 1, 100, true );

  A->AddBinScheme("pt");
  A->BinScheme("pt")->BuildBin( "Max", 0.5 ); // pT<0.5 GeV
  A->BinScheme("pt")->BuildBin( "Min", 0.5 ); // pT>0.5 GeV

  // perform the analysis ==================================
  A->Execute();

  // for reference, here is a print out of HistosDAG
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
