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
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 30000; // use this to limit the number of events
  A->writeSimpleTree = true; // write SimpleTree (for one bin)

  // set binning scheme ====================================

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Full");
  A->BinScheme("y")->BuildBin("Min",0.03);
  A->BinScheme("y")->BuildBin("Min",0.05);
  A->BinScheme("y")->BuildBin("Min",0.10);

  // final states
  A->AddFinalState("pipTrack");
  //A->AddFinalState("pimTrack");
  //A->AddFinalState("KpTrack");
  //A->AddFinalState("KmTrack");
  A->AddFinalState("jet");

  // set reconstruction method
  A->SetReconMethod("Ele");


  // perform the analysis ==================================
  A->Execute();
};
