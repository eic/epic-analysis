// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)
#include "AnalysisDelphes.h"

// cross section in Q2 bins
void analysis_xsecQ(
    TString configFile="datarec/arc/crossCheck*.root", /* delphes tree(s) */
    TString outfilePrefix="xsecQ" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      configFile,
      outfilePrefix
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

  // Q2 bins
  A->BinScheme("q2")->BuildBins(10,1,121,false);

  // y minima
  A->BinScheme("y")->BuildBin("Min",0.03);
  A->BinScheme("y")->BuildBin("Min",0.05);
  A->BinScheme("y")->BuildBin("Min",0.10);


  // perform the analysis ==================================
  A->Execute();
};
