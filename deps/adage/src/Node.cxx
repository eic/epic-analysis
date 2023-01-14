#include "Node.h"

ClassImp(Node)

using std::cout;
using std::cerr;
using std::endl;

// constructor
Node::Node(Int_t nodeType_, TString id_, CutDef *cut_, Int_t binNum_)
  : nodeType(nodeType_)
  , id(id_)
  , cut(cut_)
  , binNum(binNum_)
  , active(true)
  , debug(false)
{
  UnstageOps();
};

// print info for this node
void Node::Print() {
  TString nodeTypeStr;
  switch(nodeType) {
    case NT::bin:     nodeTypeStr="Bin";          break;
    case NT::control: nodeTypeStr="Control";      break;
    case NT::root:    nodeTypeStr="Root";         break;
    case NT::leaf:    nodeTypeStr="Leaf";         break;
    case NT::multi:   nodeTypeStr="MultiControl"; break;
    default:          nodeTypeStr="Unknown";
  };
  cout << "NODE ::: " << id << " ::: " << nodeTypeStr << endl;
  cout << "  Inputs:" << endl;
  for(Node *N : inputList) cout << "    - " << N->GetID() << endl;
  cout << "  Outputs:" << endl;
  for(Node *N : outputList) cout << "    - " << N->GetID() << endl;
  //cout << "  Temps:" << endl;
  //for(Node *N : tempList) cout << "    - " << N->GetID() << endl;
};


// add inputs and outputs to nodes 
// - inputs: nodes that connect to this node via incoming arrows
// - outputs: nodes that connect to this node via outgoing arrows
// - attempts to add a duplicate input or output will only print a warning, unless `silence` is true
void Node::AddInput(Node *N, Bool_t silence) { AddIO(N,inputList,"input",silence); };
void Node::AddOutput(Node *N, Bool_t silence) { AddIO(N,outputList,"output",silence); };
void Node::AddIO(Node *N, std::vector<Node*> &list, TString listName, Bool_t silence) {
  if(debug) {
    cout << "ADD TO NODE: " << listName << " " << N->GetID() << " to " << this->GetID() << endl;
  };
  for(Node *M : list) {
    if( M->GetID() == N->GetID() ) {
      if(!silence) cerr << "WARNING: tried to add duplicate node "
                        << M->GetID() << " to " << listName << " list" << endl;
      return;
    };
  };
  list.push_back(N);
};


// remove inputs or outputs
void Node::RemoveInput(Node *N) { inputList.erase(std::remove(inputList.begin(),inputList.end(),N)); };
void Node::RemoveOutput(Node *N) { outputList.erase(std::remove(outputList.begin(),outputList.end(),N)); };


// reset a control node by unstaging lambdas
void Node::UnstageOps() {
  this->StageInboundOp([](Node *N,NodePath *P){});
  this->StageOutboundOp([](Node *N,NodePath *P){});
};

// other
TString Node::GetCutType() { return cut->GetCutType(); };
TString Node::GetVarName() { 
  switch(GetNodeType()) {
    case NT::bin:
      if(cut) return cut->GetVarName(); else return "";
      break;
    case NT::control: 
    case NT::multi:
      return outputList.at(0)->GetVarName(); // return VarName of the bins this node controls
      break;
    default: return "";
  };
};

// conditional controls
// - control for inbound lambda:
void Node::ConditionalControl(bool B) {
  if(B) return; // do nothing if B==true
  if(GetNodeType()!=NT::control && GetNodeType()!=NT::multi) return;
  // disconnect all outputs, and store them in `tempList`
  tempList.clear();
  for(Node *N : outputList) tempList.push_back(N);
  outputList.clear();
};
// - control for outbound lambda
void Node::EndConditionalControl() {
  // reconnect all outputs, from `tempList`
  for(Node *N : tempList) AddOutput(N,true);
  tempList.clear();
};


Node::~Node() {
  if(cut) delete cut;
  auto del_vec = [] (auto v) {
    for(auto e : v) if(e) delete e;
    v.clear();
  };
  del_vec(inputList);
  del_vec(outputList);
  del_vec(tempList);
};
