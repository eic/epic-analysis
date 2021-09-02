R__LOAD_LIBRARY(Largex)
#include "Analysis.h"

/* template analysis macro
 * - runs an analysis and nothing more (e.g., it does not define bins)
 * - useful if you want to use your own analysis code
 */
void analysis_template(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.template" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 10000; // use this to limit the number of events


  // set binning scheme ====================================
  /* do nothing -> single bin */


  // perform the analysis ==================================
  A->Execute();
};
