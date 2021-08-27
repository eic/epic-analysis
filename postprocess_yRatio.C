R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// plot ratios of distributions for various y minima
void postprocess_yRatio(
    TString infile="out/yRatio.dire_5x41.brian.hiDiv.root"
) {

  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  // - specify settings such as diagonalized binnings
  Analysis *A = new Analysis();

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // find which bin has full-y range ============================
  Int_t byFull = -1;
  for(int by : P->GetBinNums("y")) {
    if(P->GetBinCut("y",by)->GetCutType()=="Full") {
      byFull = by;
      break;
    };
  };
  if(byFull<0) {
    cerr << "ERROR: no full y-cut bin found" << endl;
    return;
  };


  // loop over bins ==================================================
  // TODO: when Analysis::histSet is generalized, loops like this will
  // be significantly simpler
  for(int bpt : P->GetBinNums("pt")) {
  for(int bx  : P->GetBinNums("x")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int bq  : P->GetBinNums("q2")) {
  for(int bfs : P->GetBinNums("finalState")) {

    // loop over other y bins, drawing ratios of each of them to the full-y bin
    for(int by : P->GetBinNums("y")) {
      if(by==byFull) continue;

      // ALGORITHM: draw ratios of two Histos
      P->DrawRatios(
          Form("yRatio%d",by),
          A->GetHistosName(bpt,bx,bz,bq,by,bfs),
          A->GetHistosName(bpt,bx,bz,bq,byFull,bfs)
          );
    };

    // finish ALGORITHM - called after the loop to save canvases
    P->FinishDrawRatios("yRatio_"+A->GetHistosName(bpt,bx,bz,bq,byFull,bfs));

  }}}}};

  // finish
  cout << P->GetOutfileName() << " written" << endl;
  P->Finish();
};
