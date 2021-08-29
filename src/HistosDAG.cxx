#include "HistosDAG.h"

ClassImp(HistosDAG)

HistosDAG::HistosDAG()
  : debug(true)
{
  InitializeDAG();
};


// payload operator, executed on the specified Histos object
void HistosDAG::ForEach(std::function<void(Histos*)> op) {
  // TODO: implement GetHistos
  //auto mainOp = [op,GetHistos](NodePath P){ op(GetHistos(P)); };
  auto testOp = [](NodePath P){
    std::cout << "ALGORITHM on ";
    Node::PrintPath(P);
  };
  Payload(testOp);
};


HistosDAG::~HistosDAG() {
};

