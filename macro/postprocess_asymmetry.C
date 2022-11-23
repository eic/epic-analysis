// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Duane Byer, Christopher Dilks

R__LOAD_LIBRARY(Sidis-eic)

void postprocess_asymmetry(
    TString infile="out/asymmetry.dire_5x41.brian.hiDiv.root"
) {

  // setup
  PostProcessor *P = new PostProcessor(infile);

  // payload operator: execute PostProcessor::DrawSingle for each bin defined in analysis_asymmetry.C
  P->Op()->Payload( [&P](Histos *H) {
    P->DrawSingle(H,"phiH","");
    P->DrawSingle(H,"phiS","");
    P->DrawSingle(H,"phiSivers","");
    P->DrawSingle(H,"phiCollins","");
  });

  // execution and cleanup
  P->Execute();
  P->Finish();
};
