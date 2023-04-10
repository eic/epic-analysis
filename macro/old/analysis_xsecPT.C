// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)
#include "AnalysisDelphes.h"

// cross section in pT bins
void analysis_xsecPT(
    TString configFile="datarec/arc/crossCheck*.root", /* delphes tree(s) */
    TString outfilePrefix="xsecPT" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      configFile,
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
