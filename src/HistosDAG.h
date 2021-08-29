#ifndef HistosDAG_
#define HistosDAG_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"

// largex-eic
#include "DAG.h"
#include "Histos.h"


class HistosDAG : public DAG
{
  public:
    HistosDAG();
    ~HistosDAG();

    // payload operator, executed on the specified Histos object
    void ForEach(std::function<void(Histos*)> op);

  private:
    Bool_t debug;

  ClassDefOverride(HistosDAG,1);
};

#endif
