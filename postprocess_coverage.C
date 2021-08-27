R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_coverage(
    TString infile="out/coverage.dire_5x41.brian.hiDiv.root"
) {

  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  // - specify settings such as diagonalized binnings
  Analysis *A = new Analysis();

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // loop over bins ==================================================
  // TODO: when Analysis::histSet is generalized, loops like this will
  // be significantly simpler
  for(int bpt : P->GetBinNums("pt")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int by : P->GetBinNums("y")) {
  for(int bfs : P->GetBinNums("finalState")) {

    // loop over (x,Q2) grid
    for(int bx  : P->GetBinNums("x")) {
    for(int bq  : P->GetBinNums("q2")) {

      // ALGORITHM: draw a specific plot for this (x,Q2) bin
      P->DrawSingle(
          Form("coverage_q%d_x%d",bq,bx),
          A->GetHistosName(bpt,bx,bz,bq,by,bfs),
          "etaVsP"
          );
    };

    // finish ALGORITHM - called after the loop to save summary canvases
    //P->FinishGrid(); // TODO: if you want to make a grid of plots, you could
                       // implement the print out here

  }}}}};

  // finish
  cout << P->GetOutfileName() << " written" << endl;
  P->Finish();
};
