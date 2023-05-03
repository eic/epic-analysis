#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <functional>

// ROOT
#include <TSystem.h>
#include <TObject.h>
#include <TNamed.h>
#include <TString.h>

// adage
#include "CutDef.h"
#include "NodePath.h"

class NodePath;

// Node Types
namespace NT {
  enum nodeType_enum {
    bin,     /* labels a bin, contains a CutDef */
    control, /* a node that can contain lambdas */
    root,    /* top-level control node, with only outgoing edges */
    leaf,    /* bottom-level control node, with only incoming edges */
    multi,   /* a special control node that contains a lambda to alter the leaf node lambda */
  };
};


class Node : public TObject
{
  public:

    // constructor: nodeType (see NT above), and unique ID string;
    // a CutDef and bin number can be stored too, useful for bin nodes
    Node(Int_t nodeType_=NT::bin, TString id_="0", CutDef *cut_=nullptr, Int_t binNum_=-1);
    ~Node();

    // accessors and modifiers
    Int_t GetNodeType() { return nodeType; };
    void SetNodeType(Int_t nodeType_) { nodeType=nodeType_; };
    TString GetID() { return id; };
    void SetID(TString id_) { id=id_; };
    Int_t GetBinNum() { return binNum; };
    void SetBinNum(Int_t binNum_) { binNum=binNum_; };
    CutDef *GetCut() { return cut; };
    TString GetCutType();
    TString GetVarName();
    void Print(); // print info for this node


    // inputs and outputs accessors and modifiers
    // - inputs: nodes that connect to this node via incoming arrows
    // - outputs: nodes that connect to this node via outgoing arrows
    void AddInput(Node *N, Bool_t silence=false);
    void AddOutput(Node *N, Bool_t silence=false);
    std::vector<Node*> GetInputs() { return inputList; };
    std::vector<Node*> GetOutputs() { return outputList; };
    void RemoveInput(Node *N);
    void RemoveOutput(Node *N);

    // active state: use this to attach a boolean to a node; for example,
    // you can enable and disable certain nodes prior to a traversal
    void SetActiveState(Bool_t active_) { active=active_; };
    Bool_t IsActive() { return active; };

    // conditional controls: allow/disallow DAG depth-first traversal to pass
    // through a control node
    // - first call `ConditionalControl(bool B)` anywhere in the inbound
    //   lambda; if `B==true` then nothing happens, but if `B==false`, all of
    //   this control node's outputs will be removed and stored in `tempList`,
    //   which forces the depth-first traversal to call the outbound lambda
    //   immediately after the inbound lambda (see `DAG::ExecuteOps`); thus 
    //   we replicate the following behavior:
    //   ```
    //       for(Outerloop) {
    //         InboundOperator()
    //         if(B) {
    //           for(Subloop) ...
    //         }
    //         OutboundOperator()
    //       }
    //       ```
    // - If `B==false`, the disconnected edges will not be automatically
    //   reconnected; therefore you likely want to call `EndConditionalControl()`
    //   in the outbound lambda, if `ConditionalControl` was called in the
    //   inbound lambda; `EndConditionalControl()` will re-connect the disconnected
    //   nodes
    // - if you choose not to call `EndConditionalControl()`, any disconnection
    //   remains until you do (which could be useful in some cases)
    // - it may not be useful to call `ConditionalControl()` in the outbound
    //   lambda, since its effect will only occur the next time the traversal
    //   tries to descend through the node
    void ConditionalControl(bool B);
    void EndConditionalControl();


    // lambda operators: attach lamdas to a node (usu. control nodes),
    // and provide the ability to execute them
    // - inbound lambda is executed when traversal enters the node
    template<class O> void StageInboundOp(O op) { inboundOp = FormatOp(op); };
    // - outbound lambda is executed when traversal is ready to backtrack
    template<class O> void StageOutboundOp(O op) { outboundOp = FormatOp(op); };
    // - clear both inbound and outbound lambdas
    void UnstageOps();
    // - lambda execution
    void ExecuteInboundOp(NodePath *P) { inboundOp(this,P); };
    void ExecuteOutboundOp(NodePath *P) { outboundOp(this,P); };


    // format lambdas with proper arguments, to allow easy overloading
    static std::function<void(Node*,NodePath*)> FormatOp(std::function<void(Node*,NodePath*)> op) {
      return op;
    };
    static std::function<void(Node*,NodePath*)> FormatOp(std::function<void(NodePath*,Node*)> op) {
      return [op](Node *N, NodePath *P){ op(P,N); };
    };
    static std::function<void(Node*,NodePath*)> FormatOp(std::function<void(NodePath*)> op) {
      return [op](Node *N, NodePath *P){ op(P); };
    };
    static std::function<void(Node*,NodePath*)> FormatOp(std::function<void(Node*)> op) {
      return [op](Node *N, NodePath *P){ op(N); };
    };
    static std::function<void(Node*,NodePath*)> FormatOp(std::function<void()> op) {
      return [op](Node *N, NodePath *P){ op(); };
    };


  protected:
    void AddIO(Node *N, std::vector<Node*> &list, TString listName, Bool_t silence);


  private:
    CutDef *cut;
    Bool_t debug;
    Int_t nodeType;
    TString id;
    Int_t binNum;
    Bool_t active;
    std::vector<Node*> inputList;
    std::vector<Node*> outputList;
    std::vector<Node*> tempList;
    std::function<void(Node*,NodePath*)> inboundOp;
    std::function<void(Node*,NodePath*)> outboundOp;


  ClassDef(Node,1);
};
