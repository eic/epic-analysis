// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

/* draw Q2 distribution, useful to check Q2 weights
 */
void q2_distribution(
    TString configFile="datarec/delphes/5x41/delphes.config",   TString outfilePrefix="q2dist.delphes.5x41"
    // TString configFile="datarec/delphes/10x100/delphes.config", TString outfilePrefix="q2dist.delphes.10x100"
    // TString configFile="datarec/delphes/18x275/delphes.config", TString outfilePrefix="q2dist.delphes.18x275"
    // TString configFile="datarec/ecce/22.1/5x41/files.config",   TString outfilePrefix="q2dist.ecce.5x41"
    // TString configFile="datarec/ecce/22.1/10x100/files.config", TString outfilePrefix="q2dist.ecce.10x100"
    // TString configFile="datarec/ecce/22.1/18x275/files.config", TString outfilePrefix="q2dist.ecce.18x275"
) {

  // setup analysis
  Analysis *A;
  if (configFile.Contains("athena"))
    A = new AnalysisAthena(configFile, outfilePrefix);
  else if(configFile.Contains("ecce"))
    A = new AnalysisEcce(configFile, outfilePrefix);
  else
    A = new AnalysisDelphes(configFile, outfilePrefix);

  // settings
  A->SetReconMethod("Ele");
  // A->maxEvents = 10000; // limiter

  // cuts
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // final states
  A->AddFinalState("pipTrack");

  // analyze
  A->Execute();

  // draw Q2 distributions, in both linear and log scale
  auto P = new PostProcessor("out/"+outfilePrefix+".root");
  auto Draw = [&P] (Histos *H) {
    auto hist = H->Hist("Q2");
    auto canv = new TCanvas("q2canv","q2canv",1600,800);
    canv->Divide(2,1);
    for(int pad=1; pad<=2; pad++) {
      canv->cd(pad);
      canv->GetPad(pad)->SetGrid(1,1);
      canv->GetPad(pad)->SetLogx();
      if(pad==2) canv->GetPad(pad)->SetLogy();
      hist->Draw();
    }
    canv->SaveAs(P->GetPngDir()+"/q2dist.png");
  };
  P->Op()->Payload(Draw);
  P->Execute();
  P->Finish();
}
