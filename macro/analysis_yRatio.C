R__LOAD_LIBRARY(Largex)

// ratios of histograms with y-cut enabled to those with y-cut disabled
void analysis_yRatio(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="yRatio" /* output filename prefix*/
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
  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state
  //A->AddFinalState("jet"); // jets

  // set binning scheme ====================================
  // z ranges
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
  A->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Full"); // a bin with no y-cut
  A->BinScheme("y")->BuildBin("Min",0.03);
  A->BinScheme("y")->BuildBin("Min",0.05);
  A->BinScheme("y")->BuildBin("Min",0.10);

  // perform the analysis ==================================
  A->Execute();
};
