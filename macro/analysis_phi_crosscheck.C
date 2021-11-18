R__LOAD_LIBRARY(Largex)

void analysis_phi_crosscheck(
    TString infiles="datarec/tutorial.config", /* list of input files */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=-25, /* crossing angle [mrad] */
    TString outfilePrefix="phi.crosscheck" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state

  // define cuts ====================================
  // common cuts:
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)
  // tighter cuts:
  A->AddBinScheme("eta"); A->BinScheme("eta")->BuildBin("Range",-1.0,1.0); // central
  //A->AddBinScheme("x"); A->BinScheme("x")->BuildBin("Min",0.05);

  // set binning scheme ====================================
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBins( 7, 0.2, 0.9 );

  // perform the analysis ==================================
  A->Execute();
};
