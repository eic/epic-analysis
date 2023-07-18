// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

R__LOAD_LIBRARY(EpicAnalysis)

void analysis_dihadron(
    TString configFile="datarec/dihadron.test/10x100/files.config", 
    TString outfilePrefix="dihadron.test" /* output filename prefix*/
) {

  //outfilePrefix+="_DA";
  // setup analysis ========================================
  AnalysisEpic *A = new AnalysisEpic(
      configFile,
      outfilePrefix
      );
  
  //  A->maxEvents = 1000; // use this to limit the number of events
  A->SetReconMethod("ele"); // set reconstruction method
    
  A->AddFinalState("pippimDihadron"); // two-pion dihadron
  
  A->includeOutputSet["2h"] = true; // Dihadron final state variables

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("dihadron_xF1"); A->BinScheme("dihadron_xF1")->BuildBin("Min",0.0); // xF1 > 0 (first hadron is in CFR)
  A->AddBinScheme("dihadron_xF2"); A->BinScheme("dihadron_xF2")->BuildBin("Min",0.0); // xF2 > 0 (second hadron is in CFR)
  A->AddBinScheme("dihadron_pTlab1");  A->BinScheme("dihadron_pTlab1")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit for first hadron)
  A->AddBinScheme("dihadron_pTlab2");  A->BinScheme("dihadron_pTlab2")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit for first hadron)

  // perform the analysis ==================================
  A->Execute();

};
