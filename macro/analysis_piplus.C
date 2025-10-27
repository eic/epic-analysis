// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

R__LOAD_LIBRARY(EpicAnalysis)

void analysis_piplus(
    TString configFile="datarec/epic.25.08.0/18x275/files.config", /* delphes tree(s) */
    TString outfilePrefix="piplus.18x275" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisEpic *A = new AnalysisEpic(
      configFile,
      outfilePrefix
      );

  //  A->maxEvents = 500; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  A->writeSidisTree = true;

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.05,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.05,0.95); // 0.05 < z < 0.95
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)



  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
