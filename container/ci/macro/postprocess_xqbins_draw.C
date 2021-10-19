R__LOAD_LIBRARY(Largex)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_xqbins_draw(
    TString infile="out/macro.test.root"
) {

  PostProcessor *P = new PostProcessor(infile);

  P->Op()->Payload( [&P](Histos *H) {
      P->DrawSingle(H,"Q2vsX","COLZ");
      P->DrawSingle(H,"etaVsP","COLZ");
      });

  P->Execute();
  P->Finish(); 
};
