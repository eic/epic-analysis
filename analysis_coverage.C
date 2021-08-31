R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

// eta vs. p in bins of (x,Q2)
void analysis_coverage(
    TString infiles="datarec/pythia8NCDISS3_10x100_Q21_cross-0.025.root", /* delphes tree(s) */
    Double_t eleBeamEn=10, /* electron beam energy [GeV] */
    Double_t ionBeamEn=100, /* ion beam energy [GeV] */
    Double_t crossingAngle=-0.025, /* crossing angle [mrad] */
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
  A->BinScheme("q2")->BuildBins( 5, 1,    1000,  true);
  A->BinScheme("x")->BuildBins ( 10, 1e-4, 1,    true );

  A->BinScheme("eta")->BuildBin("Range", -1, -2);
  A->BinScheme("eta")->BuildBin("Range", -2, -3.5);
  A->BinScheme("eta")->BuildBin("Range", 1, 2);
  A->BinScheme("eta")->BuildBin("Range", 2, 3.5);

  A->BinScheme("eta")->BuildBin("Range", -1, 1);
  

  A->BinScheme("p")->BuildBin("Max", 2.0);
  A->BinScheme("p")->BuildBin("Max", 4.0);
  A->BinScheme("p")->BuildBin("Max", 10.0);
  A->BinScheme("p")->BuildBin("Min", 10.0);


  
  // perform the analysis ==================================
  A->Execute();
};
