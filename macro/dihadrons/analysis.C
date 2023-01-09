// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

void analysis(
    TString source="delphes",
    Double_t eleBeamEn=5,  Double_t ionBeamEn=41
    // Double_t eleBeamEn=18, Double_t ionBeamEn=275
    )
{
  Analysis *A;
  TString configFile, outfilePrefix;
  if(source=="delphes") {
    configFile = Form("datarec/delphes/%dx%d/delphes.config",(int)eleBeamEn,(int)ionBeamEn);
    outfilePrefix = Form("dihadrons.delphes.%dx%d",(int)eleBeamEn,(int)ionBeamEn);
    A = new AnalysisDelphes(configFile, outfilePrefix);
  } else {
    fmt::print(stderr,"ERROR: source '{}' not implemented\n",source);
    return;
  }

  A->maxEvents = 30000; // limiter
  A->writeSimpleTree = true;
  A->SetReconMethod("Ele");
  A->includeOutputSet["1h"] = true; // optionally output single-hadron plots

  // dihadron final state ==================================
  A->includeOutputSet["2h"] = true; // include the output set
  A->AddFinalState("pipTrack");     // call `AddFinalState` exactly 2 times: once for each hadron
  A->AddFinalState("pimTrack");

  // cuts ==================================================
  // - inclusive cuts
  A->AddBinScheme("w");     A->BinScheme("w")->BuildBin("Min",3.0);         // W > 3 GeV
  A->AddBinScheme("y");     A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  // - single-hadron cuts:
  /* don't use these for the dihadron final state, since they will only apply to the second hadron
  A->AddBinScheme("z");     A->BinScheme("z")->BuildBin("Range",0.2,0.9);   // 0.2 < z < 0.9
  A->AddBinScheme("xF");    A->BinScheme("xF")->BuildBin("Min",0.0);        // xF > 0
  A->AddBinScheme("ptLab"); A->BinScheme("ptLab")->BuildBin("Min",0.1);     // pT_lab > 0.1 GeV (tracking limit)
  */

  A->Execute();
};
