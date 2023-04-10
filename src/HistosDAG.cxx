// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

#include "HistosDAG.h"

ClassImp(HistosDAG)

// build the DAG from specified bin scheme
void HistosDAG::Build(std::map<TString,BinSet*> binSchemes) {

  // build the DAG, given the bin scheme
  // - require `"finalState"` bin to be the first layer
  BuildDAG(binSchemes,{"finalState"});

  // leaf operator, to create Histos objects
  LeafOp([this](NodePath *P){
    TString histosN = CreatePayloadName(P);
    TString histosT = CreatePayloadTitle("",P);
    if(debug) {
      std::cout << "At path " << P->PathString() << ": ";
      std::cout << "Create " << histosN << std::endl;
    };
    // instantiate Histos object
    Histos *H = new Histos(histosN,histosT);
    // add CutDefs to Histos object
    for(Node *N : P->GetBinNodes()) { H->AddCutDef(N->GetCut()); };
    // append to `histosMap`
    InsertPayloadData(P,H);
  });

  // execution
  if(debug) std::cout << "Begin Histos instantiation..." << std::endl;
  ExecuteAndClearOps();
};
