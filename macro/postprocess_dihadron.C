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
        P->DrawSingle(H,"dihadron_phiH","HIST");
        P->DrawSingle(H,"dihadron_phiRperp","HIST");
        P->DrawSingle(H,"dihadron_phiRT","HIST");
        P->DrawSingle(H,"dihadron_theta","HIST");
        P->DrawSingle(H,"dihadron_z","HIST");
        P->DrawSingle(H,"dihadron_z1","HIST");
        P->DrawSingle(H,"dihadron_z2","HIST");
        P->DrawSingle(H,"dihadron_xF","HIST");
        P->DrawSingle(H,"dihadron_xF1","HIST");
        P->DrawSingle(H,"dihadron_xF2","HIST");
        P->DrawSingle(H,"dihadron_Mh","HIST");
        P->DrawSingle(H,"dihadron_Mx","HIST");
        P->DrawSingle(H,"dihadron_pLab","HIST");
        P->DrawSingle(H,"dihadron_pTlab","HIST");
        P->DrawSingle(H,"dihadron_pT","HIST");
      }
  );

  P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  P->Execute();
  
  // finish ===================================================
  P->Finish(); 

};
