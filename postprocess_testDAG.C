R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// test DAG implementation
void postprocess_testDAG(
    TString infile="out/yRatioDAG.dire_5x41.brian.hiDiv.root"
) {
  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // operators ====================================================

  auto testOp = [](NodePath *P){
    cout << "BIN DATA:" << endl;
    for(Node *N : P->GetBinNodes()) {
      cout << " - " << N->GetID()
        << ":  " << N->GetCut()->GetCutTitle()
        << endl;
    };
  };

  P->Op()->BeforeSubloop(
    //{"x"},
    //{"q2"},
    {"x","q2"},
    //{"q2","x"},
    //{"x","q2","y"},
    //{"Q2"}, // typo
    //{"x,q2"}, // typo
    //[](NodePath *p){ cout << "NODEPATH: " << p->PathString() << endl; }
    //[](NodePath *p){ cout << "BINLIST: " << p->BinListString() << endl; }
    //[](NodePath *p){ cout << "Cuts: " << p->CutListString() << endl; }
    testOp
  );

  /*
  P->Op()->AfterSubloop({"q2"},[](NodePath *p){
    cout << "BINLIST: " << p->BinListString() << endl;
  });
  */

  P->Op()->PrintBreadth("DAG breadth traversal:");
    
  P->Op()->ForEach([&P](Histos* H){
    P->DrawSingle(H,"Q2vsX","COLZ");
  });

  // execution ====================================================
  P->Execute();

  // finish =====================================================
  cout << P->GetOutfileName() << " written" << endl;
  P->Finish();
};
