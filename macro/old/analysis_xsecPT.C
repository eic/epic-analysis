R__LOAD_LIBRARY(Sidis-eic)
#include "AnalysisDelphes.h"

// cross section in pT bins
void analysis_xsecPT(
    TString infiles="datarec/arc/crossCheck*.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="xsecPT" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );


  // set binning scheme ====================================

  // bin 1
  A->BinScheme("x")->BuildBin("CenterDelta", 0.3, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.7, 0.05 );
  A->BinScheme("q2")->BuildBin("Range", 4.0, 9.0 );

  // bin 2
  A->BinScheme("x")->BuildBin("CenterDelta", 0.1, 0.05 );
  A->BinScheme("z")->BuildBin("CenterDelta", 0.7, 0.05 );
  A->BinScheme("q2")->BuildBin("Range", 4.0, 9.0 );

  // diagonalization of (x,z,Q2) bins
  A->diagonalXZQ = true;

  // pT bins
  A->BinScheme("pt")->BuildBins(10,1e-2,3,true);

  // y minima
  A->BinScheme("y")->BuildBin("Min",0.03);
  A->BinScheme("y")->BuildBin("Min",0.05);
  A->BinScheme("y")->BuildBin("Min",0.10);


  // perform the analysis ==================================
  A->Execute();
};
