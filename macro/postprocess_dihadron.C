// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

R__LOAD_LIBRARY(EpicAnalysis)

// make dihadron kinematics plots
void postprocess_dihadron(
    TString infile="out/dihadron.test.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // lambdas =====================================================

  // payload: draw a few plots, using PostProcessor::DrawSingle
  P->Op()->Payload(
      [&P](Histos *H) {
        P->DrawSingle(H,"dihadron_phiH","");
      }
  );

  P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  P->Execute();
  
  // finish ===================================================
  P->Finish(); 

};
