R__LOAD_LIBRARY(Largex)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_dd4hep_draw(
    TString infile="out/fullsim.test.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // lambdas =====================================================
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

  // execution ===================================================
  P->Execute();
  
  // finish ===================================================
  P->Finish(); 

};
