R__LOAD_LIBRARY(Largex)

// make resolution plots
void analysis_barak(
    TString infiles="datarec/dis-18x275-xm25.config", /* list of input files */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=25, /* crossing angle [mrad] */
    TString outfilePrefix="barak_dis-18x275-xm25" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );
  A->NBINS = 1; // use this to set the number of bins along each axis, e.g., z binning (except resolution axes) for each overall bin in e.g. x and Q2
  A->NBINSRES = 100; // use this to set the number of bins along the resolution axes for each overall bin in e.g. x and Q2
  // A->RESHIGH = 1;
  // A->RESLOW = -1;
  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  A->AddBinScheme("q2"); A->BinScheme("q2")->BuildBins( 25, 0.1,    10000,  true );
  A->AddBinScheme("x");  A->BinScheme("x")->BuildBins(  25, 0.00001, 1,    true );

  // perform the analysis ==================================
  A->Execute();
};
