// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// test drawing some single (x,Q2) bin plots, as opposed to grids of plots in
// (x,Q2) bins
// - includes drawing Q2vsX plots, as a simple sanity check for the binning
void postprocess_bin_test(TString infile) {
  PostProcessor *P = new PostProcessor(infile);
  P->Op()->Payload( [&P](Histos *H) {
      P->DrawSingle(H,"Q2vsX","COLZ");
      P->DrawSingle(H,"etaVsP","COLZ");
      });
  P->Execute();
  P->Finish(); 
};
