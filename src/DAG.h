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

    void Initialize();

    DAGnode *GetNode(TString id_, Bool_t silence=false);

    DAGnode *GetTopNode();
    DAGnode *GetBottomNode();
    DAGnode *GetUniqueNode(Int_t type_,TString typeStr);

    void AddNode(DAGnode *N, Bool_t silence=false);
    void RenameNode(DAGnode *N, TString newName, Int_t newType=-1);

    void AddLayer(std::vector<DAGnode*> nodes);

    void AddEdge(TString inID, TString outID);
    void AddEdge(DAGnode *inN, TString outID);
    void AddEdge(TString inID, DAGnode *outN);
    void AddEdge(DAGnode *inN, DAGnode *outN);

    void TraverseLayers(DAGnode *N, void (*lambda)(DAGnode*));

    void Print(TString header="DAG");

  protected:
    Bool_t Visited(TString id_);

  private:
    Bool_t debug;
    std::map<TString,DAGnode*> nodeMap;
    std::vector<TString> visitList;

  ClassDef(DAG,1);
};

#endif
