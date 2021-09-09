R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

// eta vs. p in bins of (x,Q2)
void analysis_coverage(
    TString infiles="datarec/dire_5x41.brian.hiDiv.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="coverage" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );


  // set binning scheme ====================================

  // (x,Q2) grid of bins, equal width in log-scale
  A->BinScheme("q2")->BuildBins( 3, 1,    100,  true);
  A->BinScheme("x")->BuildBins ( 3, 5e-2, 1,    true );

  // perform the analysis ==================================
  A->Execute();
};
