// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// ratios of histograms with y-cut enabled to those with y-cut disabled
void analysis_yRatio(
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

  // define cuts ===========================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  // z ranges
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
  A->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Max",0.95); // a bin with no minimum y-cut
  A->BinScheme("y")->BuildBin("Range",0.03,0.95);
  A->BinScheme("y")->BuildBin("Range",0.05,0.95);
  A->BinScheme("y")->BuildBin("Range",0.10,0.95);

  // perform the analysis ==================================
  A->Execute();
};
