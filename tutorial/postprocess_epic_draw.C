// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_epic_draw(
    TString infile="out/tutorial.epic.root"
) {

  // setup postprocessor ========================================
  // - this will read in the Histos objects and your defined
  //   binning, to construct a HistosDAG object
  PostProcessor *P = new PostProcessor(infile);

  // print DAG ==================================================
  // - here we print out the HistosDAG object, to show the structure
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  // - observe how the structure changes for the different bins
  //   you specified in the analysis macro
  P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // lambdas =====================================================
  // - we now define what operations we want to happen while
  //   traversing through the DAG; this example shows how to define
  //   a payload operator, which will act on every multidimensional
  //   bin

  // payload: draw a few plots, using PostProcessor::DrawSingle
  // - this lambda expression will be executed on every multidimensional
  //   bin's Histos object
  // - capture the pointer to PostProcessor `P` by reference, so we
  //   have access to it while the lambda executes
  // - lambda arguments can include a Histos pointer, and optinally
  //   a `NodePath` pointer, which gives binning information; here
  //   we just need the Histos pointer
  // - `PostProcessor::Op()` is just an interface (alias) for the
  //   `HistosDAG` object; see classes `HistosDAG` and its parent `DAG`
  //   for available methods
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

  // print DAG ============================
  // - we only added a payload, we did not restructure the DAG,
  //   therefore this second printout should have the same structure
  //   as the first
  P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  // - after you have defined all your operators, run them:
  P->Execute();
  
  // finish ===================================================
  // - PostProcessor::Finish() is NECESSARY to close output streams
  // - output file names and directories will be printed to stdout
  P->Finish(); 

};
