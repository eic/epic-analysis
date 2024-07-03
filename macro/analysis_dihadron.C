// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

R__LOAD_LIBRARY(EpicAnalysis)


struct WeightsTest : public WeightsSivers {
  Double_t Asymmetry(Double_t x, Double_t z, Double_t Q2, Double_t pt) const override {
    return 0.068;
  }
};



void analysis_dihadron(
    TString configFile="datarec/6.26.2024.test/10x100/files.config", 
    //TString configFile="datarec/latest.6.11.2024___epic.24.05.0_5x41/5x41/files.config", 
    //TString configFile="datarec/single_piplus_study_14_per_Q2___epic.23.10.0_10x100/10x100/files.config",
    TString outfilePrefix="dihadron.6.26.2024.nov2023" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisEpic *A = new AnalysisEpic(
      configFile,
      outfilePrefix
      );

//   Weights* weights = new WeightsSum({
//       new WeightsUniform(),
//       new WeightsTest()
//     });
//   A->SetDihadronWeights(weights);   

  //  A->maxEvents = 1000; // use this to limit the number of events
  A->SetReconMethod("ele"); // set reconstruction method
  A->writeDiSidisTree = true; // save Dihadron kinematics to simple TTree

  A->AddFinalState("pippimDihadron"); // two-pion dihadron
  
  A->includeOutputSet["1h"] = false; // Single hadron final state variables
  A->includeOutputSet["2h"] = true; // Dihadron final state variables

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.05,0.95); // 0.05 < y < 0.95
  A->AddBinScheme("dihadron_xF1"); A->BinScheme("dihadron_xF1")->BuildBin("Min",0.0); // xF1 > 0 (first hadron is in CFR)
  A->AddBinScheme("dihadron_xF2"); A->BinScheme("dihadron_xF2")->BuildBin("Min",0.0); // xF2 > 0 (second hadron is in CFR)
  A->AddBinScheme("dihadron_pTlab1");  A->BinScheme("dihadron_pTlab1")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit for first hadron)
  A->AddBinScheme("dihadron_pTlab2");  A->BinScheme("dihadron_pTlab2")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit for first hadron)

  // perform the analysis ==================================
  A->Execute();

};
