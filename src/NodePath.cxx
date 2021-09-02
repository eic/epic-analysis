#include "NodePath.h"

ClassImp(NodePath)

// constructor
NodePath::NodePath()
{
};

// return a string of bin names
TString NodePath::BinListString() {
  TString ret = "[";
  for(auto it : GetBinNodes()) ret += " " + it->GetID();
  ret += " ]";
  return ret;
};

// print a string of node names
void NodePath::PrintBinList() { std::cout << BinListString() << std::endl; };


// get string listing of node names
TString NodePath::PathString() {
  TString ret = "[";
  for(auto it : nodes) ret += " " + it->GetID();
  ret += " ]";
  return ret;
};

// print string listing of node names
void NodePath::PrintPath() { std::cout << PathString() << std::endl; };

// return set of bin nodes only (mainly for HistosDAG::GetHistos())
std::set<Node*> NodePath::GetBinNodes() {
  std::set<Node*> ret;
  for(auto N : nodes) {
    if(N->GetNodeType()==NT::bin) ret.insert(N);
  };
  return ret;
};

// return list of cuts
TString NodePath::CutListString() {
  TString ret = "";
  Bool_t first = true;
  for(auto N : GetBinNodes()) {
    if(!first) ret += ",  ";
    first=false;
    ret += N->GetCut()->GetCutTitle();
  };
  return ret;
};




NodePath::~NodePath() {
};
