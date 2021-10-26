R__LOAD_LIBRARY(Largex)

void postprocess_phi_crosscheck(
    TString infile="out/phi.crosscheck.root"
) {

  // cleanup old image files
  gROOT->ProcessLine(".! rm -v out/phi.crosscheck.images/*.png");

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // lambdas =====================================================
  P->Op()->Payload(
      [&P](Histos *H) {
        P->DrawSingle(H,"phiH","hist");
      });

  //P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  P->Execute();
  P->Finish(); 
};
