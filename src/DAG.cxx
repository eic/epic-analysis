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


void DAG::Initialize() {
  nodeMap.clear();
  AddEdge(new DAGnode(tTop,"top0"),new DAGnode(tBottom,"bottom0"));
};


void DAG::AddNode(DAGnode *N, Bool_t silence) {
  TString id_ = N->GetID();
  if(debug) cout << "ADD TO DAG: Node " << id_ << endl;
  if(GetNode(id_,true)) {
    if(!silence)
      cerr << "WARNING: tried to add duplicate node " << id_ << " to DAG" << endl;
    return;
  } else {
    nodeMap.insert(std::pair<TString,DAGnode*>(id_,N));
  };
};


void DAG::RenameNode(DAGnode *N, TString newName, Int_t newType) {
  nodeMap.erase(N->GetID());
  N->SetID(newName);
  if(newType>=0) N->SetNodeType(newType);
  nodeMap.insert(std::pair<TString,DAGnode*>(newName,N));
};


void DAG::AddLayer(std::vector<DAGnode*> nodes) {
  auto C = GetBottomNode();
  RenameNode(C,"tmp",tControl);
  AddNode(new DAGnode(tBottom,"bottom0"));
  for(auto N : nodes) {
    AddEdge(C,N);
    AddEdge(N,"bottom0");
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


void DAG::Print(TString header) {
  TString sep = "=========================";
  cout << endl << header << endl << sep << endl;
  TraverseLayers( GetTopNode(), [](DAGnode *N){ N->Print(); cout << endl; });
  cout << sep << endl;
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


DAGnode *DAG::GetTopNode() { return GetUniqueNode(tTop,"top"); };
DAGnode *DAG::GetBottomNode() { return GetUniqueNode(tBottom,"bottom"); };
DAGnode *DAG::GetUniqueNode(Int_t type_,TString typeStr) {
  DAGnode * ret = nullptr;
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

