R__LOAD_LIBRARY(Sidis-eic)

// analysis in bins of (p,eta)
void analysis_p_eta(
    TString infiles="datarec/canyonlands-v1.2/5x41/files.config", // default, for manual local testing
    Double_t eleBeamEn=5,
    Double_t ionBeamEn=41,
    Double_t crossingAngle=-25,
    TString outfilePrefix="coverage.fullsim",
    TString reconMethod="Ele"
) {

  // setup analysis ========================================
  Analysis *A;
  if(outfilePrefix.Contains("fullsim"))
       A = new AnalysisAthena(  infiles, eleBeamEn, ionBeamEn, crossingAngle, outfilePrefix );
#ifndef EXCLUDE_DELPHES
  else A = new AnalysisDelphes( infiles, eleBeamEn, ionBeamEn, crossingAngle, outfilePrefix );
#endif

  A->SetReconMethod(reconMethod); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  // A->writeSimpleTree = true;

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
