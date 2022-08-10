#include "NodePath.h"

ClassImp(NodePath)

// constructor
NodePath::NodePath()
{
};

// return a string of bin names -- for pretty print
TString NodePath::BinListString() {
  TString ret = "[";
  for(auto it : GetSortedBins()) ret += " " + it->GetID();
  ret += " ]";
  return ret;
};

// return a string of bin names -- for object names
TString NodePath::BinListName() {
  TString ret = "";
  for(auto it : GetSortedBins()) ret += "_" + it->GetID();
  ret(TRegexp("^_")) = "";
  return ret;
};

// print a string of node names
void NodePath::PrintBinList() { std::cout << BinListString() << std::endl; };


// get string listing of node names
TString NodePath::PathString() {
  TString ret = "[";
  for(auto it : GetSortedPath()) ret += " " + it->GetID();
  ret += " ]";
  return ret;
};

// print string listing of node names
void NodePath::PrintPath() { std::cout << PathString() << std::endl; };

// return set of bin nodes only
std::set<Node*> NodePath::GetBinNodes() {
  std::set<Node*> ret;
  for(auto N : nodes) { if(N->GetNodeType()==NT::bin) ret.insert(N); };
  return ret;
};

// return list of bin nodes in DAG order
std::vector<Node*> NodePath::GetSortedBins() {
  std::vector<Node*> ret;
  for(auto N : GetSortedPath()) { if(N->GetNodeType()==NT::bin) ret.push_back(N); };
  return ret;
};


// return list of cuts
TString NodePath::CutListString() {
  TString ret = "";
  Bool_t first = true;
  for(auto N : GetSortedBins()) {
    if(!first) ret += ",  ";
    first=false;
    ret += N->GetCut()->GetCutTitle();
  };
  return ret;
};


// return the node that points to `N`
Node *NodePath::GetPreviousNode(Node *N, bool silence) {
  for(auto inN : N->GetInputs()) { if(nodes.find(inN)!=nodes.end()) return inN; };
  if(!silence) std::cerr << "ERROR: node " << N->GetID() << " has no previous node in this path" << std::endl;
  return nullptr;
};
// return the node that `N` points to
Node *NodePath::GetNextNode(Node *N, bool silence) {
  for(auto outN : N->GetOutputs()) { if(nodes.find(outN)!=nodes.end()) return outN; };
  if(!silence) std::cerr << "ERROR: node " << N->GetID() << " has no next node in this path" << std::endl;
  return nullptr;
};
// return the first or last node in the path
Node *NodePath::GetFirstNode() {
  for(auto N : nodes) { if(GetPreviousNode(N,true)==nullptr) return N; };
  return nullptr;
};
Node *NodePath::GetLastNode() {
  for(auto N : nodes) { if(GetNextNode(N,true)==nullptr) return N; };
  return nullptr;
};


// returns path in DAG order (NodePath encapsulates an unordered set of nodes)
std::vector<Node*> NodePath::GetSortedPath() {
  std::vector<Node*> ret;
  Node *N = GetFirstNode();
  while(N) {
    ret.push_back(N);
    N = GetNextNode(N,true);
  };
  return ret;
};

// return the bin node for the given variable name
Node *NodePath::GetBinNode(TString varName) {
  for(auto N : GetBinNodes()) { if(varName == N->GetVarName()) return N; };
  std::cerr << "ERROR: bin node for variable " << varName << " not found in path" << std::endl;
  return nullptr;
};

NodePath::~NodePath() {
};
