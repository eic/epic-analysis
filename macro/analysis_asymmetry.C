R__LOAD_LIBRARY(Largex)

// test new SimpleTree class
void analysis_asymmetry(
    TString infiles="datarec/arc/main_subset_10/*.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="asym" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 10000;
  A->writeSimpleTree = true;

  // set binning scheme ====================================
  /* do nothing -> single bin */


  // perform the analysis ==================================
  A->Execute();
};
