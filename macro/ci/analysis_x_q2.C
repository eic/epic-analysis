R__LOAD_LIBRARY(Largex)

// analysis in bins of (x,Q2)
void analysis_x_q2(
    TString infiles="datarec/canyonlands-v1.2/5x41/files.config", // default, for manual local testing
    Double_t eleBeamEn=5,
    Double_t ionBeamEn=41,
    Double_t crossingAngle=-25,
    TString outfilePrefix="coverage.fullsim"
) {

  // setup analysis ========================================
  Analysis *A;
  if(outfilePrefix.Contains("fullsim"))
       A = new AnalysisDD4hep(  infiles, eleBeamEn, ionBeamEn, crossingAngle, outfilePrefix );
  else A = new AnalysisDelphes( infiles, eleBeamEn, ionBeamEn, crossingAngle, outfilePrefix );

  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state

  // define cuts ===========================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  Int_t nx, nq;
  if(outfilePrefix.Contains("bin-test")) { nx=3; nq=3; } else { nx=6; nq=4; };
  A->AddBinScheme("x");  A->BinScheme("x")->BuildBins(  nx, 0.001, 1,    true );
  A->AddBinScheme("q2"); A->BinScheme("q2")->BuildBins( nq, 1,     3000,  true );

  // perform the analysis ==================================
  A->Execute();
};
