#include "DAG.h"

ClassImp(DAG)

using std::cout;
using std::cerr;
using std::endl;

// constructor
DAG::DAG()
  : debug(true)
{
};


void DAG::AddNode(DAGnode *N, Bool_t silence) {
  TString id_ = N->GetID();
  if(debug) cout << "ADD TO DAG: Node " << id_ << endl;
  if(GetNode(id_,true)) {
    if(!silence)
      cerr << "WARNING: tried to add duplicate node " << id_ << "to DAG" << endl;
    return;
  } else {
    nodeMap.insert(std::pair<TString,DAGnode*>(id_,N));
  };
};


void DAG::AddEdge(TString inID, TString outID) { this->AddEdge(this->GetNode(inID),this->GetNode(outID)); };
void DAG::AddEdge(DAGnode *inN, TString outID) { this->AddEdge(inN,this->GetNode(outID)); };
void DAG::AddEdge(TString inID, DAGnode *outN) { this->AddEdge(this->GetNode(inID),outN); };
void DAG::AddEdge(DAGnode *inN, DAGnode *outN) {
  if(inN && outN) {
    inN->AddOutput(outN);
    outN->AddInput(inN);
    this->AddNode(inN,true);
    this->AddNode(outN,true);
  } else {
    cerr << "ERROR: tried to a edge to non-existing node(s)" << endl;
  };
};


void DAG::Print() {
  cout << "\nDAG\n====================\n";
  TraverseLayers( GetTopNode(), [](DAGnode *N){ N->Print(); cout << endl; });
  cout << "====================\n";
};


DAGnode *DAG::GetNode(TString id_, Bool_t silence) {
  auto kv = nodeMap.find(id_);
  if(kv != nodeMap.end()) return kv->second;
  else {
    if(!silence)
      cerr << "ERROR: cannot find node " << id_ << " in DAG" << endl;
    return nullptr;
  };
};


DAGnode *DAG::GetTopNode() {
  DAGnode * ret = nullptr;
  for(auto kv : nodeMap) {
    auto N = kv.second;
    if(N->GetNodeType()==tTop) {
      if(ret!=nullptr) cerr << "WARNING: this DAG has more than one top node" << endl;
      ret = N;
    };
  };
  if(ret==nullptr) cerr << "WARNING: this DAG has no top node" << endl;
  return ret;
};



void DAG::TraverseLayers(DAGnode *N, void (*lambda)(DAGnode*)) {
  TString id_ = N->GetID();
  if(N->GetNodeType()==tTop) {
    lambda(N);
    visitList.clear();
  };
  if(!Visited(id_)) {
    visitList.push_back(id_);
    for(auto M : N->GetOutputs()) {
      if(!Visited(M->GetID())) lambda(M);
    };
    for(auto M : N->GetOutputs()) TraverseLayers(M,lambda);
  };
};

Bool_t DAG::Visited(TString id_) {
  return std::find(visitList.begin(),visitList.end(),id_) != visitList.end();
};


DAG::~DAG() {
};

