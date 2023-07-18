// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

R__LOAD_LIBRARY(EpicAnalysis)

void analysis_dihadron(
    TString configFile="datarec/dihadron.test/18x275/files.config", 
    TString outfilePrefix="dihadron.test" /* output filename prefix*/
) {

  //outfilePrefix+="_DA";
  // setup analysis ========================================
  AnalysisEpic *A = new AnalysisEpic(
      configFile,
      outfilePrefix
      );
  
  A->maxEvents = 1000; // use this to limit the number of events
  A->SetReconMethod("ele"); // set reconstruction method
  
  
  //  A->AddFinalState("pipTrack"); // pion final state
  //  A->AddFinalState("pimTrack"); // pion final state
  A->AddFinalState("pippimDihadron"); // two-pion dihadron
  
  A->includeOutputSet["2h"] = true; // Dihadron final state variables

  // define cuts ====================================
  //  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  //  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  //  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================
  // z ranges
  //A->AddBinScheme("z");
  //  A->BinScheme("z")->BuildBin("Min",0.2); // needed?

  // y minima
  //  A->AddBinScheme("y");
  //  A->BinScheme("y")->BuildBin("Full"); // a bin with no y-cut

  // perform the analysis ==================================
  A->Execute();

};
