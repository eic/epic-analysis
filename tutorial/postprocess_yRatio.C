// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// plot ratios of distributions for various y minima
void postprocess_yRatio(
    TString infile="out/yRatio.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // print DAG ==================================================
  P->Op()->PrintBreadth("HistosDAG Initial Setup");


  // lambdas =====================================================

  // we will need to define two payloads; this will be handled
  // by the `MultiPayload` staging function below

  // payload 1: find the full-y bin; set Hfull to that Histos pointer
  // - `Hfull` is captured by reference, since this lambda sets it
  // - Histos pointer `H` is included as an argument, since this
  //   is a payload operator
  // - the NodePath must also be included as an argument, so we 
  //   can check if the cut type of the 1D y-bin is "Full"
  Histos *Hfull = nullptr;
  auto findFullBin = [&Hfull](Histos *H, NodePath *bins) {
    if(bins->GetBinNode("y")->GetCutType() == "Max") Hfull = H;
  };

  // payload 2: draw ratio of each plot with a set y-minimum to that with
  // no y-minimum (full-y bin)
  // - `Hfull` is captured by reference, since the previous lambda
  //   determined it
  // - PostProcessor pointer `P` is also captured
  // - we only need the Histos pointer in the argument
  auto drawRatios = [&Hfull,&P](Histos *H) {
    // make sure we have a denominator (full-y bin)
    if(Hfull==nullptr) {
      cerr << "ERROR: no full-y bin found" << endl;
      return;
    };
    // skip the full-y bin (don't bother making ratio full / full )
    if(H==Hfull) return;
    // PostProcessor::DrawRatios will take the ratio of each 1D
    // histogram in the specified Histos objects
    P->DrawRatios(
        "ratio_" + Hfull->GetSetName(), // output name
        H, // numerator Histos
        Hfull // denominator Histos
        );
  };

  // after subloop over y bins
  // - finish ratio drawing, for saving canvases etc.
  auto drawRatiosFinish = [&Hfull,&P]() {
    P->FinishDrawRatios("yRatio_"+Hfull->GetSetName());
  };


  // lambda staging ==========================================
  // run MultiPayload on a subloop over y-bins, for the two lambdas
  P->Op()->MultiPayload(
      {"y"}, // list of bin layers in subloop (std::vector<TString>)
      findFullBin // payload operator
      ); 
  P->Op()->MultiPayload(
      {"y"}, // list of bin layers
      drawRatios, // payload operator
      [](){}, // empty before operator
      drawRatiosFinish // after operator
      );
  
  // note: do not define `MultiPayload` for different subloops; the 
  // reason is because `MultiPayload` creates a control node which
  // has an inbound lambda that overwrites the current payload. 
  // `MutliPayload` will also overwrite any payload operator you
  // have staged with `Payload`. You can either define one `Payload`
  // or several `MultiPayloads` on a particular subloop.


  // print DAG ============================
  // - compare the configured DAG with the initial DAG; notice
  //   the new multi-control nodes before the y bins
  P->Op()->PrintBreadth("HistosDAG Final Setup");
  //P->Op()->PrintDepth(); // print depth-first node paths

  // execution ====================================================
  P->Execute();

  // finish =====================================================
  P->Finish();
};
