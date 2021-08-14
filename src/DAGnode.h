#ifndef DAGnode_
#define DAGnode_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"


// largex-eic
#include "CutDef.h"

enum nodeType_enum {
  tBin,
  tControl,
  tTop,
  tBottom
};

class DAGnode : public TObject
{
  public:

    DAGnode(Int_t nodeType_=tBin, TString id_="0");
    ~DAGnode();

    Int_t GetNodeType() { return nodeType; };
    void SetNodeType(Int_t nodeType_) { nodeType=nodeType_; };
    TString GetID() { return id; };

    void AddInput(DAGnode *N);
    void AddOutput(DAGnode *N);

    void Print();

  protected:
    void AddNodeToList(DAGnode *N, std::vector<DAGnode*> &list, TString listName);

  private:
    Bool_t debug;
    Int_t nodeType;
    TString id;
    std::vector<DAGnode*> inputList;
    std::vector<DAGnode*> outputList;

  ClassDef(DAGnode,1);
};

#endif
