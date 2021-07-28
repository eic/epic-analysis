// test macro, to try out the new Analysis class
R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

void analysis_test() {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      "datarec/cross_5x41_25_newcard.root",
      5,
      41,
      25
      );


  // set binning scheme ====================================
  A->diagonalBinsOnly = true;
  
  // slide 11
  A->BinScheme("x")->BuildBin("CenterDelta", 0.3, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.7, 0.05 );
  A->BinScheme("pt")->BuildBin("CenterDelta", 0.5, 0.05 );

  // slide 13
  A->BinScheme("x")->BuildBin("CenterDelta", 0.6, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.5, 0.05 );
  A->BinScheme("pt")->BuildBin("CenterDelta", 0.55, 0.05 );
  // slide 14
  A->BinScheme("x")->BuildBin("CenterDelta", 0.1, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.7, 0.05 );
  A->BinScheme("pt")->BuildBin("CenterDelta", 0.15, 0.05 );

  // Q bins
  A->BinScheme("q")->BuildBins(10,1,11,false);

  // y minima
  A->BinScheme("y")->BuildBin("Min",0.05);


  // jet test, minimum pT cut
  A->BinScheme("pt_jet")->BuildBin("Min",10.0);

  // perform the analysis ==================================
  A->Execute();
};
