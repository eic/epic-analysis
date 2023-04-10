// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

void postprocess_phi_crosscheck(
    TString infile="out/phi.crosscheck.root"
) {

  // cleanup old image files
  gROOT->ProcessLine(".! rm -v out/phi.crosscheck.images/*.png");

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // lambdas =====================================================
  P->Op()->Payload(
      [&P](Histos *H) {
        P->DrawSingle(H,"phiH","hist");
      });

  //P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  P->Execute();
  P->Finish(); 
};
