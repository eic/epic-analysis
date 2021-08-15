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
// - bin: labels a bin, contains a CutDef
// - control: a node that can contain lambdas
// - root: top-level control node, with only outgoing edges
// - leaf: bottom-level control node, with only incoming edges
namespace NT {
  enum nodeType_enum {
    bin,
    control,
    root,
    leaf
  };
};


class Node : public TObject
{
  public:

    // constructor: nodeType (see NT above), and unique ID string
    Node(Int_t nodeType_=NT::bin, TString id_="0");
    ~Node();

    // accessors and modifiers
    Int_t GetNodeType() { return nodeType; };
    void SetNodeType(Int_t nodeType_) { nodeType=nodeType_; };
    TString GetID() { return id; };
    void SetID(TString id_) { id=id_; };

    // inputs and outputs accessors and modifiers
    // - inputs: nodes that connect to this node via incoming arrows
    // - outputs: nodes that connect to this node via outgoing arrows
    void AddInput(Node *N);
    void AddOutput(Node *N);
    std::vector<Node*> GetInputs() { return inputList; };
    std::vector<Node*> GetOutputs() { return outputList; };
    void RemoveInput(Node *N);
    void RemoveOutput(Node *N);

    // print info for this node
    void Print();


  protected:
    void AddIO(Node *N, std::vector<Node*> &list, TString listName);

  private:
    Bool_t debug;
    Int_t nodeType;
    TString id;
    std::vector<Node*> inputList;
    std::vector<Node*> outputList;

  ClassDef(Node,1);
};

#endif
