R__LOAD_LIBRARY(Sidis-eic)

void analysis_depolarization(
    // Double_t eleBeamEn=5,
    // Double_t ionBeamEn=41,
    /**/
    Double_t eleBeamEn=18,
    Double_t ionBeamEn=275,
    /**/
    Double_t crossingAngle=-25
    )
{

  TString infiles = Form("datarec/delphes.%dx%d.config",(int)eleBeamEn,(int)ionBeamEn);
  TString outfilePrefix = Form("depol.%dx%d",(int)eleBeamEn,(int)ionBeamEn);

  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state

  /// SIDIS common cuts
  A->AddBinScheme("w");     A->BinScheme("w")->BuildBin("Min",3.0);         // W > 3 GeV
  A->AddBinScheme("y");     A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");     A->BinScheme("z")->BuildBin("Range",0.2,0.9);   // 0.2 < z < 0.9
  A->AddBinScheme("xF");    A->BinScheme("xF")->BuildBin("Min",0.0);        // xF > 0
  A->AddBinScheme("ptLab"); A->BinScheme("ptLab")->BuildBin("Min",0.1);     // pT_lab > 0.1 GeV (tracking limit)

  /// additional cuts
  A->AddBinScheme("x"); A->BinScheme("x")->BuildBin("CenterDelta", 0.3, 0.1 );

  /// run
  A->Execute();
};
