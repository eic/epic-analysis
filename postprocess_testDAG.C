R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// test DAG implementation
void postprocess_testDAG(
    TString infile="out/yRatioDAG.dire_5x41.brian.hiDiv.root"
) {
  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // operators ====================================================

  // `controlOp` is a lambda which prints some information about the 
  // subloop control node; it is meant to demonstrate some actions
  // you can take at that point in the bin loop
  auto controlOp = [](NodePath *bins, Node *controlNode, TString vars){
    cout << endl;
    cout << "CONTROL " << vars << endl;
    cout << "  NODE = " << controlNode->GetID() << endl;
    cout << "  PATH = " << bins->PathString() << endl;
    cout << "  PREVIOUS = " << bins->GetPreviousNode(controlNode)->GetID()
         << endl;
  };

  // examples of useful subloop operators ///////////////
  auto testOp = [](NodePath *P){
    cout << endl;
    cout << "BIN DATA:" << endl;
    for(Node *N : P->GetBinNodes()) {
      cout << " - " << N->GetID()
        << ":  " << N->GetCut()->GetCutTitle()
        << endl;
    };
  };


  // stage operators on DAG ========================================

  /*
  P->Op()->BeforeSubloop(
    //{"x"},
    {"q2"},
    //{"x","q2"},
    //{"q2","x"},
    //{"x","q2","y"},
    //{"Q2"}, // typo
    //{"x,q2"}, // typo
    //[](NodePath *bins){cout << "NODES: " << bins->PathString() << endl; }
    //[](NodePath *bins){cout << "BINS: " << bins->BinListString() << endl; }
    //[](NodePath *bins){cout << "Cuts: " << bins->CutListString() << endl; }
    testOp
  );
  P->Op()->AfterSubloop({"q2"},[](NodePath *bins){
    cout << "BINLIST: " << bins->BinListString() << endl << endl;
  });
  */

  // subloops in series /////////////////////////
  // - both include the same kinematics, causing a control node to appear
  //   before the `x` layer
  // - the second `Subloop` operator will overwrite the first, which may
  //   not be useful to do; control nodes in series are not allowed,
  //   since it is possible to just combine the lambdas
  // - consider instead using `MultiPayload` if you want separate subloops
  /*
  P->Op()->Subloop({"x","q2"},[](NodePath *bins){ cout << "\nPATH one: " << bins->PathString() << endl; });
  P->Op()->Subloop({"x","q2"},[](NodePath *bins){ cout << "\nPATH two: " << bins->PathString() << endl; });
  */

  // nested subloops //////////////////////////////
  // - the first listed bin variable becomes the "control" node; all other
  //   listed variables will be included in the subloop, in order listed
  // - you can declare smaller subloops after bigger subloops, similarly
  //   to nested for loops; you can do it the other way around, but you
  //   may overwrite some subloop actions
  // - you can also make subloops over each variable, one at a time, and
  //   they will be iterated over in order
  // - as shown above, duplicate control nodes overwrite
  // - turn on/off the various subloops, to gain intuition how they behave
  // --- example 1: separate subloops
  //P->Op()->Subloop( {"y"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"y"); });
  //P->Op()->Subloop( {"x","q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x,q2"); });
  // --- example 2: nest loop over q2 inside loop over (x,q2)
  //P->Op()->Subloop( {"x","q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x,q2"); });
  //P->Op()->Subloop( {"q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"q2"); });
  // --- example 3: notice the looping is equivalent to example 2
  //P->Op()->Subloop( {"x"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x"); });
  //P->Op()->Subloop( {"q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"q2"); });
  // --- sandbox: turn on/off various subloops, to test behavior
  //P->Op()->Subloop( {"y","x","q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"y,x,q2"); });
  //P->Op()->Subloop( {"x","q2","y"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x,q2,y"); });
  //P->Op()->Subloop( {"x","q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x,q2"); });
  //P->Op()->Subloop( {"q2","x"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"q2,x"); });
  //P->Op()->Subloop( {"q2"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"q2"); });
  //P->Op()->Subloop( {"x"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"x"); });
  //P->Op()->Subloop( {"y"}, [&controlOp](NodePath *B,Node *C){ controlOp(B,C,"y"); });

  // multiple payloads //////////////////////////////
  // - if you would like to have more than one payload operator, use `MultiPayload`
  // - the first argument is a list of layers, which will define the subloop; this subloop
  //   will be run multiple times
  // - you can define as many `MultiPayload` payloads as you want
  // - if you define one `MultiPayload`, it will completely override any `Payload`
  // - do not define `MultiPayload` for more than one subloop; only one subloop will run the
  //   multiple payloads
  // - you can also include `Before` and `After` operators
  /*
  P->Op()->MultiPayload( {"x","q2"}, [](Histos *H){ cout << "MULTI-PAYLOAD 1: " << H->GetSetName() << endl; });
  P->Op()->MultiPayload( {"x","q2"}, [](Histos *H){ cout << "MULTI-PAYLOAD 2: " << H->GetSetName() << endl; });
  P->Op()->MultiPayload( {"x","q2"}, [](Histos *H){ cout << "MULTI-PAYLOAD 3: " << H->GetSetName() << endl; },
      [](){ cout < "BEFORE MULTI-PAYLOAD 3" << endl; },
      [](){ cout < "AFTER MULTI-PAYLOAD 3" << endl; }
      );
  */


  // payload
  P->Op()->Payload([&P](Histos* H){
    cout << "PAYLOAD: " << H->GetSetName() << endl;
    //P->DrawSingle(H,"Q2vsX","COLZ");
  });

  P->Op()->PrintBreadth("DAG breadth traversal:");
  P->Op()->PrintDepth("DAG depth traversal:");

  // execution ====================================================
  P->Execute();

  // finish =====================================================
  cout << endl;
  P->Finish();
};
