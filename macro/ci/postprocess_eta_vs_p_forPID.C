// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// plot eta vs. p in one bin, for PID coverage requests
void postprocess_eta_vs_p_forPID(TString infile="out/coverage.fastsim.root") {
  PostProcessor *P = new PostProcessor(infile);
  P->Op()->Payload( [&P](Histos *H) {
      P->DrawSingle(H,"etaVsP","COLZ");
      });
  P->Execute();
  P->Finish(); 
}
