#ifndef Node_
#define Node_

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

// Node Types
namespace NT {
  enum nodeType_enum {
    bin,     /* labels a bin, contains a CutDef */
    control, /* a node that can contain lambdas */
    root,    /* top-level control node, with only outgoing edges */
    leaf     /* bottom-level control node, with only incoming edges */
  };
};


class Node : public TObject
{
  public:

    // constructor: nodeType (see NT above), and unique ID string;
    // a CutDef can be stored too, useful for bin nodes
    Node(Int_t nodeType_=NT::bin, TString id_="0", CutDef *cut_=nullptr);
    ~Node();

    // accessors and modifiers
    Int_t GetNodeType() { return nodeType; };
    void SetNodeType(Int_t nodeType_) { nodeType=nodeType_; };
    TString GetID() { return id; };
    void SetID(TString id_) { id=id_; };
    CutDef *GetCut() { return cut; };
    TString GetVarName() { return ( cut==nullptr ? "" : cut->GetVarName() ); };

    // inputs and outputs accessors and modifiers
    // - inputs: nodes that connect to this node via incoming arrows
    // - outputs: nodes that connect to this node via outgoing arrows
    void AddInput(Node *N, Bool_t silence=false);
    void AddOutput(Node *N, Bool_t silence=false);
    std::vector<Node*> GetInputs() { return inputList; };
    std::vector<Node*> GetOutputs() { return outputList; };
    void RemoveInput(Node *N);
    void RemoveOutput(Node *N);

    // print info for this node
    void Print();


  protected:
    void AddIO(Node *N, std::vector<Node*> &list, TString listName, Bool_t silence);

  private:
    CutDef *cut;
    Bool_t debug;
    Int_t nodeType;
    TString id;
    std::vector<Node*> inputList;
    std::vector<Node*> outputList;

  ClassDef(Node,1);
};

#endif
