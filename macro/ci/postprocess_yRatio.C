// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// plot ratios of distributions for various y minima
void postprocess_yRatio(TString infile) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // print DAG ==================================================
  P->Op()->PrintBreadth("HistosDAG Initial Setup");


  // lambdas =====================================================

  // payload 1: find the full-y bin; set Hfull to that Histos pointer
  Histos *Hfull = nullptr;
  auto findFullBin = [&Hfull](Histos *H, NodePath *bins) {
    if(bins->GetBinNode("y")->GetCutType() == "Max") Hfull = H;
  };

  // payload 2: draw ratio of each plot with a set y-minimum to that with
  // no y-minimum (full-y bin)
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
  
  P->Op()->PrintBreadth("HistosDAG Final Setup");
  //P->Op()->PrintDepth(); // print depth-first node paths

  // execution ====================================================
  P->Execute();

  // finish =====================================================
  P->Finish();
};
