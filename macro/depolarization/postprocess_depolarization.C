// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

R__LOAD_LIBRARY(Sidis-eic)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_depolarization(
    TString infile="out/depol.delphes.5x41.root"
    // TString infile="out/depol.delphes.18x275.root"
) {

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
      });

  P->Execute();
  P->Finish(); 

};
