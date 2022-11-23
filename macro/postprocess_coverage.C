// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

R__LOAD_LIBRARY(Sidis-eic)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_coverage(
    TString infile="out/coverage.example_5x41.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // lambdas =====================================================

  // payload: draw a few plots, using PostProcessor::DrawSingle
  P->Op()->Payload(
      [&P](Histos *H) {
        P->DrawSingle(H,"Q2vsX","COLZ");
        P->DrawSingle(H,"etaVsP","COLZ");
        //P->DrawSingle(H,"x_Res","");
        //P->DrawSingle(H,"x_RvG","COLZ");
        //P->DrawSingle(H,"phiH_RvG","COLZ");
        //P->DrawSingle(H,"phiS_RvG","COLZ");
      }
      );

  P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  P->Execute();
  
  // finish ===================================================
  P->Finish(); 

};
