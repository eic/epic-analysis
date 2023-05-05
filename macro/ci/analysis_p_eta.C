// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// analysis in bins of (p,eta)
void analysis_p_eta(
    TString configFile,
    TString outfilePrefix,
    TString reconMethod="Ele"
    )
{

  // setup analysis ========================================
  Analysis *A;
  if   (outfilePrefix.Contains("epic"))     A = new AnalysisEpic(   configFile, outfilePrefix );
  else if(outfilePrefix.Contains("athena")) A = new AnalysisAthena( configFile, outfilePrefix );
  else if(outfilePrefix.Contains("ecce"))   A = new AnalysisEcce(   configFile, outfilePrefix );
#ifndef EXCLUDE_DELPHES
  else A = new AnalysisDelphes( configFile, outfilePrefix );
#endif

  A->SetReconMethod(reconMethod); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  // A->writeSidisTree = true;

  // define cuts ===========================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  A->AddBinScheme("p");   A->BinScheme("p")->BuildBins( 6, 0.1, 100, true );
  A->AddBinScheme("eta"); A->BinScheme("eta")->BuildBins( 4, -4, 4, false );

  // perform the analysis ==================================
  A->Execute();
};
