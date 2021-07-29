R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

// cross section in Q bins
void analysis_xsecQ() {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      "datarec/arc/dire_5x41.brian.hiDiv.root", /* delphes tree */
      5, /* electron beam energy [GeV] */
      41, /* ion beam energy [GeV] */
      0 /* crossing angle [mrad] */
      );


  // set binning scheme ====================================

  // bin 1 (cf Harut slide 11)
  A->BinScheme("x")->BuildBin("CenterDelta", 0.3, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.7, 0.05 );
  A->BinScheme("pt")->BuildBin("CenterDelta", 0.5, 0.05 );

  // bin 2
  A->BinScheme("x")->BuildBin("CenterDelta", 0.3, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.4, 0.05 );
  A->BinScheme("pt")->BuildBin("CenterDelta", 0.15, 0.05 );

  // diagonalization of (pt,x,z) bins
  A->diagonalPtXZ = true;

  // Q bins
  A->BinScheme("q")->BuildBins(10,1,11,false);

  // y minima
  A->BinScheme("y")->BuildBin("Min",0.03);
  A->BinScheme("y")->BuildBin("Min",0.05);
  A->BinScheme("y")->BuildBin("Min",0.10);


  // perform the analysis ==================================
  A->Execute();
};
