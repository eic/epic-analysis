R__LOAD_LIBRARY(Largex)

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
