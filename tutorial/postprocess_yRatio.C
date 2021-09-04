R__LOAD_LIBRARY(Largex)

// plot ratios of distributions for various y minima
void postprocess_yRatio(
    TString infile="out/yRatio.example_5x41.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);


  // lambdas =====================================================

  // payload: find the full-y bin; set Hfull to that Histos pointer
  Histos *Hfull = nullptr;
  auto findFullBin = [&Hfull](Histos *H, NodePath *bins) {
    if(bins->GetBinNode("y")->GetCutType() == "Full") Hfull = H;
  };

  // payload: draw ratio
  auto drawRatios = [&Hfull,&P](Histos *H, NodePath *bins) {
    // make sure we have a denominator (full-y bin)
    if(Hfull==nullptr) {
      cerr << "ERROR: no full-y bin found" << endl;
      return;
    };
    // skip the full-y bin
    if(H==Hfull) return;
    // PostProcessor::DrawRatios will take the ratio of each 1D
    // histogram in the specified Histos objects
    P->DrawRatios(
        "ratio_" + bins->GetBinNode("y")->GetID(), // output name
        H, // numerator Histos
        Hfull // denominator Histos
        );
  };

  // control after: finish ratio drawing, for saving canvases etc.
  auto drawRatiosFinish = [&Hfull,&P]() {
    P->FinishDrawRatios("yRatio_"+Hfull->GetSetName());
  };


  // operator staging ==========================================
  // run MultiPayload on a subloop over y-bins, for the two lambdas
  P->Op()->MultiPayload(
      {"y"},
      findFullBin
      ); 
  P->Op()->MultiPayload(
      {"y"},
      drawRatios, // payload
      [](){}, // empty before op
      drawRatiosFinish // after op
      );


  // print DAG ============================
  P->Op()->PrintBreadth();
  P->Op()->PrintDepth();

  // execution ====================================================
  P->Execute();

  // finish =====================================================
  P->Finish();



  /*
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
  */
};
