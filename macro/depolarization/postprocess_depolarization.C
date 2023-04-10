// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_depolarization(TString infile) {

  PostProcessor *P = new PostProcessor(infile);

  P->Op()->Payload(
      [&P](Histos *H) {
        for(auto histName : H->VarNameList) {
          if(histName.Contains("epsilon") || histName.Contains("depol")) {
            if(histName.Contains("Q2vsX")) 
              P->DrawSingle(H,histName,"SURF",1);
            else
              P->DrawSingle(H,histName,"COLZ",1,true);
          }
        }
        P->GetOutfile()->cd("/");
        H->Hist("Q2vsX")->Write();
        H->Hist("y")->Write();
      });

  P->Execute();
  P->Finish(); 

};
