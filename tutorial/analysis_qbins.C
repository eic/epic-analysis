R__LOAD_LIBRARY(Largex)

// run in Q2 bins, for two pT ranges
void analysis_qbins(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.xqbins" /* output filename prefix*/
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

  // set binning scheme ====================================
  A->BinScheme("q2")->BuildBins( 10, 1, 100, true );
  A->BinScheme("pt")->BuildBin( "Max", 0.5 ); // pT<0.5 GeV
  A->BinScheme("pt")->BuildBin( "Min", 0.5 ); // pT>0.5 GeV


  // perform the analysis ==================================
  A->Execute();
};
