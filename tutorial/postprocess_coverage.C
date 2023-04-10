// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_coverage(
    TString infile="out/coverage.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  //P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // lambdas =====================================================

  /* first, let's loop through the (x,Q2) bins, for only full-range
   * (eta,p) bins; in each (x,Q2) bin we draw eta vs. p
   * - define the subloop over (x,Q2) bins
   * - include `ConditionalControl(bool B)` in the inbound lambda ("before 
   *   subloop"); the subloop will only be iterated through if `B==true`
   * - you must also include `EndConditionalControl()` in the outbound
   *   lambda ("after subloop"), in order to return the DAG to normal
   * - both `ConditionalControl` and `EndConditionalControl` are functions
   *   of the current control Node, so you must include a Node pointer in the
   *   lambda arguments
   */

  // lambda to check if all bins in the current bin path (NodePath)
  // are full; this will be an inbound lambda and will call
  // `ConditionalControl`
  auto checkIfFull = [](NodePath *bins, Node *controlNode){
    // print a header
    cout << endl;
    cout << "before subloop; bins: " << bins->BinListString() << endl;
    cout << "                cuts: " << bins->CutListString() << endl;
    // initialize boolean `B`
    bool B = true;
    // loop through bins in bin path; if one of them is not "Full", set B=false;
    // be careful not to do this check on "finalState" bins, which
    // are never going to be "Full"
    for(Node *bin : bins->GetBinNodes()) {
      if(bin->GetCutType()!="Full" && bin->GetVarName()!="finalState") {
        // we also need to ignore the bins which are just cut definitions, since they are not "Full":
        if( bin->GetVarName() != "w"
         && bin->GetVarName() != "y"
         && bin->GetVarName() != "z"
         && bin->GetVarName() != "xF"
         && bin->GetVarName() != "ptLab"
        ) {
          B = false;
          break;
        };
      };
    };
    // call `ConditionalControl` on the current node, which is
    // the subloop's control node
    controlNode->ConditionalControl(B);
    // print what we will do
    if(B) cout << "-> this is a full bin, proceed into subloop:" << endl;
    else cout << "-> this is not a full bin, cancel subloop" << endl;
    // show behavior behind the scenes:
    // print the control node status, to show inputs and outputs;
    // notice that if B==false, all this node's outputs are moved to 
    // "Temps", thus this node is disconnected and depth-first DAG
    // traversal will not continue into the subloop; otherwise if B==true
    // the outputs are not changed
    /*
    cout << "  control node status:" << endl;
    controlNode->Print();
    cout << endl;
    */
  };


  // lambda to end the conditional control; this will be an outbound
  // lambda, called after the subloop
  // - the only thing this lambda does is call `EndConditionalControl()`;
  //   everything else is for printing information to demonstrate how
  //   the node connections change
  // - the NodePath pointer `bins` is only needed for printout
  auto endIfFull = [](NodePath *bins, Node *controlNode){
    // call `EndConditionalControl()` to revert this node's connections
    controlNode->EndConditionalControl();
    // show behavior behind the scenes:
    // print the control node status, to show inputs and outputs; notice
    // that any "Temps" are moved back to "Outputs", so that the 
    // full connection is restored, in preparation for the next time
    // this node is visited
    /*
    cout << "\nafter subloop; bins: " << bins->BinListString() << endl;
    cout << "  control node status:" << endl;
    controlNode->Print();
    cout << endl;
    */
  };


  // lambdas staging =====================================================

  // print a divider
  cout << "\n\n#####################################################\n";
  cout << "etaVsP in (x,Q2) subloop\n";
  cout << "#####################################################\n\n";

  // now we define the (x,Q2) subloop, staging `checkIfFull` as the inbound
  // lambda, and an outbound lambda which only calls `EndConditionalControl`
  // - the inbound lambda only is aware of what (eta,p) bin we are in,
  //   therefore only eta and p are checked for being "Full"
  P->Op()->Subloop(
      {"x","q2"}, // layers list
      checkIfFull, // inbound lambda
      endIfFull // outbound lambda
      );

  // stage payload lambda, to draw eta vs. P plots; we only care about
  // drawing non-full (x,Q2) bins
  P->Op()->Payload([&P](Histos *H, NodePath *bins){ 
    if(bins->GetBinNode("x")->GetCutType()=="Full" || bins->GetBinNode("q2")->GetCutType()=="Full") return;
    P->DrawSingle(H,"etaVsP","COLZ");
  });

  // print the DAG and execute
  //P->Op()->PrintBreadth("HistosDAG Setup: etaVsP in (x,Q2) subloop");
  P->Execute();

  // note: by default, `PostProcessor::Execute()` clears all staged lambdas
  // ( you can call `Execute(false)` if you do not want lambdas to be cleared)


  // now let's do the "opposite": make Q2 vs. x plots in bins of (eta,p):

  // print a divider
  cout << "\n\n#####################################################\n";
  cout << "Q2vsX in (eta,p) subloop\n";
  cout << "#####################################################\n\n";

  // similar to before, we want to be sure we are in a full (x,Q2) bin before
  // the subloop over (eta,p) bins; the inbound lambda only is aware of what
  // (x,Q2) bin we are in, therefore only x and Q2 are checked for being "Full"
  P->Op()->Subloop( {"eta","p"}, checkIfFull, endIfFull );
  
  // define a similar payload, omitting full (eta,p) bins
  P->Op()->Payload([&P](Histos *H, NodePath *bins){ 
    if(bins->GetBinNode("eta")->GetCutType()=="Full" || bins->GetBinNode("p")->GetCutType()=="Full") return;
    P->DrawSingle(H,"Q2vsX","COLZ");
  });

  // execution
  //P->Op()->PrintBreadth("HistosDAG Setup: Q2vsX in (eta,p) subloop");
  P->Execute();

  cout << "\n\n#####################################################\n\n";


  
  // finish ===================================================
  P->Finish(); 

};
