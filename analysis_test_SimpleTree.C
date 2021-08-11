R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

// test new SimpleTree class
void analysis_test_SimpleTree(
    TString infiles="datarec/arc/crossCheck_5x41.run8.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="test.simple.tree" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->maxEvents = 10000;
  A->writeSimpleTree = true;

  // set binning scheme ====================================
  /* do nothing -> single bin */


  // perform the analysis ==================================
  A->Execute();
};
