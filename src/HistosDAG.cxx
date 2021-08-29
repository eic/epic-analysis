#include "HistosDAG.h"

ClassImp(HistosDAG)

HistosDAG::HistosDAG()
  : debug(true)
{
  InitializeDAG();
};

HistosDAG::~HistosDAG() {
};

