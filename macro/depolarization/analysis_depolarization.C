// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

void analysis_depolarization(
    TString configFile,
    TString outfilePrefix,
    TString upstream="epic"
    )
{

  Analysis *A;
  if(upstream=="delphes") {
    A = new AnalysisDelphes(configFile, outfilePrefix);
  } else if(upstream=="epic") {
    A = new AnalysisEpic(configFile, outfilePrefix);
  } else if(upstream=="ecce") {
    A = new AnalysisEcce(configFile, outfilePrefix);
  } else if(upstream=="athena") {
    A = new AnalysisAthena(configFile, outfilePrefix);
  } else {
    fmt::print(stderr,"ERROR: unknown upstream='{}'\n",upstream);
    return 1;
  }

  // A->maxEvents = 100000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  A->includeOutputSet["depolarization"] = true; // include depolarization plots

  /// SIDIS common cuts
  A->AddBinScheme("w");     A->BinScheme("w")->BuildBin("Min",3.0);         // W > 3 GeV
  A->AddBinScheme("y");     A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");     A->BinScheme("z")->BuildBin("Range",0.2,0.9);   // 0.2 < z < 0.9
  A->AddBinScheme("xF");    A->BinScheme("xF")->BuildBin("Min",0.0);        // xF > 0
  A->AddBinScheme("ptLab"); A->BinScheme("ptLab")->BuildBin("Min",0.1);     // pT_lab > 0.1 GeV (tracking limit)

  /// run
  A->Execute();
};
