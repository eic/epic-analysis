#include "DAG.h"

ClassImp(DAG)

using std::cout;
using std::cerr;
using std::endl;

// constructor
DAG::DAG()
  : debug(true)
{
  Initialize();
};


// initialize DAG, with only the root node connected to the leaf node
void DAG::Initialize() {
  nodeMap.clear();
  AddEdge(new Node(NT::root,"root0"),new Node(NT::leaf,"leaf0"));
};


// search for a node by ID
Node *DAG::GetNode(TString id_, Bool_t silence) {
  auto kv = nodeMap.find(id_);
  if(kv != nodeMap.end()) return kv->second;
  else {
    if(!silence)
      cerr << "ERROR: cannot find node " << id_ << " in DAG" << endl;
    return nullptr;
  };
};


// get unique nodes (viz. root and leaf)
Node *DAG::GetRootNode() { return GetUniqueNode(NT::root,"root"); };
Node *DAG::GetLeafNode() { return GetUniqueNode(NT::leaf,"leaf"); };
Node *DAG::GetUniqueNode(Int_t type_,TString typeStr) {
  Node * ret = nullptr;
  for(auto kv : nodeMap) {
    auto N = kv.second;
    if(N->GetNodeType()==type_) {
      if(ret!=nullptr) cerr << "WARNING: this DAG has more than one " << typeStr << " node" << endl;
      ret = N;
    };
  };
  if(ret==nullptr) cerr << "WARNING: this DAG has no " << typeStr << " node" << endl;
  return ret;
};


// add node to DAG
void DAG::AddNode(Int_t nodeType, TString id, Bool_t silence) {
  AddNode(new Node(nodeType,id),silence);
};
void DAG::AddNode(Node *N, Bool_t silence) {
  TString id_ = N->GetID();
  if(debug) cout << "ADD TO DAG: Node " << id_ << endl;
  if(GetNode(id_,true)) {
    if(!silence)
      cerr << "WARNING: tried to add duplicate node " << id_ << " to DAG" << endl;
    return;
  } else {
    nodeMap.insert(std::pair<TString,Node*>(id_,N));
  };
};


// add edges
void DAG::AddEdge(Node *inN, Node *outN) {
  if(inN && outN) {
    inN->AddOutput(outN);
    outN->AddInput(inN);
    this->AddNode(inN,true);
    this->AddNode(outN,true);
  } else {
    cerr << "ERROR: tried to a edge to non-existing node(s)" << endl;
  };
};


// rename or repurpose a node
void DAG::ModifyNode(Node *N, TString newName, Int_t newType) {
  nodeMap.erase(N->GetID());
  N->SetID(newName);
  if(newType>=0) N->SetNodeType(newType);
  nodeMap.insert(std::pair<TString,Node*>(newName,N));
};


// remove a node
void DAG::RemoveNode(Node *N) {
  for(auto inN : N->GetInputs()) RemoveEdge(inN,N);
  for(auto outN : N->GetOutputs()) RemoveEdge(N,outN);
  nodeMap.erase(N->GetID());
};


// remove an edge
void DAG::RemoveEdge(Node *inN, Node *outN) {
  inN->RemoveOutput(outN);
  outN->RemoveInput(inN);
};


// add a layer of nodes, fully connected to the last layer of the DAG
// - primary usage is to add a layer of bins from a BinSet
void DAG::AddLayer(std::vector<Node*> nodes) {
  auto C = GetLeafNode();
  ModifyNode(C,"tmp",NT::control);
  AddNode(NT::leaf,"leaf0");
  for(auto N : nodes) {
    AddEdge(C,N);
    AddEdge(N,GetNode("leaf0"));
  };
  Simplify();
};


// print the whole DAG (breadth first)
void DAG::Print(TString header) {
  TString sep = "=========================";
  cout << endl << header << endl << sep << endl;
  TraverseBreadth( GetRootNode(), [](Node *N){ N->Print(); cout << endl; });
  cout << sep << endl;
};


// breadth-first traversal
void DAG::TraverseBreadth(Node *N, std::function<void(Node*)> lambda) {
  TString id_ = N->GetID();
  if(N->GetNodeType()==NT::root) {
    lambda(N);
    visitList.clear();
  };
  if(!Visited(id_)) {
    visitList.push_back(id_);
    for(auto M : N->GetOutputs()) {
      if(!Visited(M->GetID())) lambda(M);
    };
    for(auto M : N->GetOutputs()) TraverseBreadth(M,lambda);
  };
};


// depth-first traversal
void DAG::TraverseDepth(Node *N, std::function<void(Node*)> lambda) {
};


// traversal helper which checks if a node has been visited
Bool_t DAG::Visited(TString id_) {
  return std::find(visitList.begin(),visitList.end(),id_) != visitList.end();
};


// simplify DAG: convert all control nodes into full connections 
// between the adjacent layers; all control nodes will be removed
void DAG::Simplify() {
  auto nodeMapCopy = nodeMap;
  for(auto kv : nodeMapCopy) RemoveControl(kv.first);
};

// remove a control node, replacing it with full connection between
// adjacent layers
void DAG::RemoveControl(TString id) { RemoveControl(GetNode(id)); };
void DAG::RemoveControl(Node *N) {
  if(N->GetNodeType()!=NT::control) return;
  for(auto inN : N->GetInputs()) {
    for(auto outN : N->GetOutputs()) AddEdge(inN,outN);
  };
  RemoveNode(N);
};


DAG::~DAG() {
};

