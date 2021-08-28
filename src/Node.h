#ifndef Node_
#define Node_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <functional>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"

// largex-eic
#include "CutDef.h"

// DAG path of Nodes is handled by std::set of Node pointers
class Node;
typedef std::set<Node*> NodePath;

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

    // lambda staging: attach lamdas to a node (usu. control nodes)
    // - lambdas can be executed in depth-first DAG traversal
    // - inbound lambda is executed when traversal enters the node
    // - outbound lambda is executed when traversal is ready to backtrack
    void StageInboundOp(std::function<void(Node*,NodePath)> op) { inboundOp=op; };
    void StageOutboundOp(std::function<void(Node*,NodePath)> op) { outboundOp=op; };
    void UnstageOps(); // reset
    // lambda execution
    void ExecuteInboundOp(NodePath P) { inboundOp(this,P); };
    void ExecuteOutboundOp(NodePath P) { outboundOp(this,P); };

    // print info for this node
    void Print();

    // helper: print a specific NodePath
    static void PrintPath(NodePath P) {
      std::cout << "[";
      for(auto it : P) std::cout << " " << it->GetID();
      std::cout << "]" << std::endl;
    };


  protected:
    void AddIO(Node *N, std::vector<Node*> &list, TString listName, Bool_t silence);

  private:
    CutDef *cut;
    Bool_t debug;
    Int_t nodeType;
    TString id;
    std::vector<Node*> inputList;
    std::vector<Node*> outputList;
    std::function<void(Node*,NodePath)> inboundOp;
    std::function<void(Node*,NodePath)> outboundOp;

  ClassDef(Node,1);
};

#endif
