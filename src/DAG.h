#ifndef DAG_
#define DAG_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"


// largex-eic
#include "CutDef.h"
#include "DAGnode.h"

class DAG : public TObject
{
  public:

    DAG();
    ~DAG();

    DAGnode *GetNode(TString id_, Bool_t silence=false);

    void AddNode(DAGnode *N, Bool_t silence=false);

    void AddEdge(TString inID, TString outID);
    void AddEdge(DAGnode *inN, TString outID);
    void AddEdge(TString inID, DAGnode *outN);
    void AddEdge(DAGnode *inN, DAGnode *outN);

    void Print();


  private:
    Bool_t debug;
    std::map<TString,DAGnode*> nodeMap;

  ClassDef(DAG,1);
};

#endif
