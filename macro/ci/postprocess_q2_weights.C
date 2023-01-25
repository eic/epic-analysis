// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// plot the Q2 distribution, to check if Q2 weights are correctly applied
void postprocess_q2_weights(TString infile="out/coverage.fastsim.root") {
  PostProcessor *P = new PostProcessor(infile);
  P->Op()->Payload( [&P](Histos *H) {
      P->DrawSingle(H,"Q2");
      });
  P->Execute();
  P->Finish(); 
}
