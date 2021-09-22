R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// plot ratios of distributions for various y minima
void postprocess_asymmetry(
    TString infile="out/asymmetry.dire_5x41.brian.hiDiv.root"
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
  for(int bx  : P->GetBinNums("x")) {
  for(int by : P->GetBinNums("y")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int bq  : P->GetBinNums("q2")) {
  for(int bfs : P->GetBinNums("finalState")) {
      P->DrawSingle(
          A->GetHistosName(bpt,bx,bz,bq,by,bfs),
          "phiS"
          );
  }}}}}};

  // finish
  cout << P->GetOutfileName() << " written" << endl;
  P->Finish();
};

