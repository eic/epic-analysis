#include "DAG.h"

ClassImp(DAG)

using std::cout;
using std::cerr;
using std::endl;

// constructor
DAG::DAG()
  : debug(false)
  , activeEvent(true)
{
  InitializeDAG();
};


// initialize DAG, with only the root node connected to the leaf node
void DAG::InitializeDAG() {
  nodeMap.clear();
  layerMap.clear();
  valueMap.clear();
  valueMapExternal.clear();
  AddEdge(new Node(NT::root,"root_0"),new Node(NT::leaf,"leaf_0"));
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
    nodeMap.insert({id_,N});
  };
};


// add edges
void DAG::AddEdge(Node *inN, Node *outN, Bool_t silence) {
  if(inN && outN) {
    inN->AddOutput(outN,silence);
    outN->AddInput(inN,silence);
    this->AddNode(inN,true);
    this->AddNode(outN,true);
  } else {
    cerr << "ERROR: tried to a edge to non-existing node(s)" << endl;
  };
};


// rename or repurpose a node
void DAG::ModifyNode(Node *N, TString newName, Int_t newType) {
  if(debug) cout << "RENAME NODE: " << N->GetID();
  nodeMap.erase(N->GetID());
  N->SetID(newName);
  if(newType>=0) N->SetNodeType(newType);
  nodeMap.insert({newName,N});
  if(debug) cout << " to " << N->GetID() << endl;
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
void DAG::AddLayer(BinSet *BS) {
  std::vector<Node*> nodes;
  int binNum=0; // enumerate bin number
  TObjArrayIter next(BS->GetBinList());
  while(CutDef *cut = (CutDef*) next()) {
    nodes.push_back(new Node( NT::bin, cut->GetVarName()+Form("_%d",binNum), cut, binNum ));
    binNum++;
  };
  AddLayer(nodes);
  layerMap.insert({BS->GetVarName(),new BinSet(*BS)});
};
// - add layer of nodes
void DAG::AddLayer(std::vector<Node*> nodes) {
  // convert the leaf node to a temporary control node
  auto C = GetLeafNode();
  ModifyNode(C,"tmp",NT::control);
  // create new leaf
  AddNode(NT::leaf,"leaf_0");
  // connect new layer after temporary control node and before new leaf
  for(auto N : nodes) {
    AddEdge(C,N);
    AddEdge(N,GetNode("leaf_0"));
  };
  RepatchToFull(C);
};

// get the BinSet associated to a layer
BinSet *DAG::GetBinSet(TString varName) {
  BinSet *BS;
  try {
    BS = layerMap.at(varName);
  } catch(const std::out_of_range &ex) {
    std::cerr << "ERROR: layer \"" << varName << "\" not found in DAG" << std::endl;
    return nullptr;
  };
  return BS;
};


// print the whole DAG (breadth first)
void DAG::PrintBreadth(TString header) {
  TString sep = "=========================";
  auto printOp = [](Node *N){ N->Print(); cout << endl; };
  cout << endl << header << endl << sep << endl;
  TraverseBreadth(GetRootNode(),printOp);
  cout << sep << endl;
};
// print the whole DAG (depth first)
void DAG::PrintDepth(TString header) {
  TString sep = "=========================";
  auto printOp = [this](Node *N,NodePath *P){ P->PrintPath(); };
  cout << endl << header << endl << sep << endl;
  TraverseDepth(GetRootNode(),printOp);
  cout << sep << endl;
};
// print unique paths from the root node to the leaf node
void DAG::PrintLeafPaths(TString header) {
  TString sep = "=========================";
  auto printOp = [this](Node *N,NodePath *P){
    if(N->GetNodeType()==NT::leaf) P->PrintPath();
  };
  cout << endl << header << endl << sep << endl;
  TraverseDepth(GetRootNode(),printOp);
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
void DAG::TraverseDepth(Node *N, std::function<void(Node*,NodePath*)> lambda, NodePath P) {
  P.nodes.insert(N);
  lambda(N,&P);
  for(auto M : N->GetOutputs()) TraverseDepth(M,lambda,P);
};


// run each node's staged lambdas, while traversing depth first; if `activeNodesOnly`, all nodes
// must have `active==true` (by default all nodes are active)
void DAG::ExecuteOps(Bool_t activeNodesOnly, Node *N, NodePath P) {
  if(N==nullptr) N = GetRootNode();
  if(activeNodesOnly && N->IsActive()==false) return;
  P.nodes.insert(N);
  N->ExecuteInboundOp(&P);
  for(auto M : N->GetOutputs()) ExecuteOps(activeNodesOnly,M,P);
  N->ExecuteOutboundOp(&P);
};


// clear all staged lambdas
void DAG::ClearOps() {
  for(auto kv : nodeMap) kv.second->UnstageOps();
};


// associate a layer (bin scheme) with a value, set by a lambda
// - add the lambda to `valueMap`, iff there exists a layer for `binVarName` 
// - if `force` is true, always add the lambda to `valueMap` (maybe slower)
void DAG::SetBinSchemeValue(
    TString binVarName,
    std::function<Double_t()> binValueLambda,
    Bool_t force
    )
{
  if( force ? true : layerMap.find(binVarName) != layerMap.end() )
    valueMap.insert({binVarName,binValueLambda});
};
// similar to `SetBinSchemeValue`, but for `valueMapExternal`, for CutDefs of type `external`
void DAG::SetBinSchemeValueExternal(
    TString binVarName,
    std::function<Bool_t(Node*)> binValueLambda,
    Bool_t force
    )
{
  if( force ? true : layerMap.find(binVarName) != layerMap.end() )
    valueMapExternal.insert({binVarName,binValueLambda});
};


// check all bins, and activate nodes for which the CutDef is satisfied
void DAG::CheckBins() {
  // activate bin nodes for which the CutDef is satisfied
  auto CheckBin = [this] (Node *N) {
    if(N->GetNodeType()==NT::bin) {
      Bool_t active;
      Double_t val;
      // if it's an external cut, check the boolean returned from the `valueMapExternal` lambda
      if(N->GetCut()->IsExternal()) {
        try {
          active = valueMapExternal.at(N->GetVarName())(N);
        } catch(const std::out_of_range &ex) {
          cerr << "\nERROR: variable \"" << N->GetVarName() << "\" not found in DAG::valueMapExternal" << endl;
          active = false;
        };
      }
      // otherwise, check the value returned from the `valueMap` lambda using `CutDef::CheckCut`
      else {
        try {
          // get value associated to this variable, and check cut
          val = valueMap.at(N->GetVarName())();
          active = N->GetCut()->CheckCut(val);
        } catch(const std::out_of_range &ex) {
          cerr << "\nERROR: variable \"" << N->GetVarName() << "\" not found in DAG::valueMap" << endl;
          active = false;
        };
      };
      N->SetActiveState(active);
    };
  };
  TraverseBreadth(CheckBin);
  // set `activeEvent` if there is at least one multidimensional bin to fill
  // - calls `ExecuteOps` with `activeNodesOnly=true`; the `LeafOp` will only
  //   execute if the depth-first traversal reaches the leaf node, which
  //   can only happen if there is at least one NodePath with all of its Nodes active
  activeEvent = false;
  LeafOp( [this](){ activeEvent = true; } );
  ExecuteOps(true);
}


// set all nodes to active (the default) or inactive (if active_==false)
void DAG::ActivateAllNodes(Bool_t active_) {
  for(auto kv : nodeMap) kv.second->SetActiveState(active_);
};


// traversal helper which checks if a node has been visited
Bool_t DAG::Visited(TString id_) {
  return std::find(visitList.begin(),visitList.end(),id_) != visitList.end();
};


// convert control patches to full patches; control node(s) will be removed
void DAG::RepatchToFull(TString id) { RepatchToFull(GetNode(id)); };
void DAG::RepatchToFull(Node *N) {
  if(N->GetNodeType()!=NT::control && N->GetNodeType()!=NT::multi) return;
  // connect each input to each output
  for(auto inN : N->GetInputs()) {
    for(auto outN : N->GetOutputs()) AddEdge(inN,outN,true);
  };
  // remove control node and its connections
  RemoveNode(N);
};
void DAG::RepatchAllToFull() {
  auto nodeMapCopy = nodeMap;
  for(auto kv : nodeMapCopy) RepatchToFull(kv.first);
};


// re-patch a layer to the current leaf node's output, and re-patch the
// adjacent layers together; the current leaf node will become a control
// node, and a new leaf node is created
void DAG::RepatchToLeaf(TString varName) {
  // check if the variable exists in this DAG
  Bool_t nodeExists = false;
  for(auto kv : nodeMap) {
    auto N = kv.second;
    if(N->GetNodeType()==NT::bin) {
      if(N->GetVarName()==varName) {
        nodeExists=true;
        break;
      };
    };
  };
  if(!nodeExists) {
    cerr << "\nERROR: variable \"" << varName << "\" not found in DAG; do you have a typo?\n\n";
    return;
  };
  // check if a control node for this variable exists; if so, repatch it to full
  TString controlID = varName+"_control";
  auto Ctmp = GetNode(controlID,true);
  if(Ctmp!=nullptr) RepatchToFull(Ctmp);
  // convert the leaf node to a control node
  auto C = GetLeafNode();
  ModifyNode(C,controlID,NT::control);
  // create new leaf
  AddNode(NT::leaf,"leaf_0");
  auto L = GetLeafNode();
  // repatching lambda
  auto repatch = [&C,&L,&varName,this](Node *N) {
    // only act on nodes in the specified layer
    if(N->GetNodeType()!=NT::bin || N->GetVarName()!=varName) return;
    // re-patch adjacent layers together
    for(auto inN : N->GetInputs()) {
      for(auto outN : N->GetOutputs()) this->AddEdge(inN,outN,true);
    };
    // disconnect node N from adjacent layers
    for(auto inN : N->GetInputs()) this->RemoveEdge(inN,N);
    for(auto outN : N->GetOutputs()) this->RemoveEdge(N,outN);
    // re-patch node N to be between control and leaf nodes
    this->AddEdge(C,N,true);
    this->AddEdge(N,L,true);
    // mark node as visited, so we don't try to repatch it again when
    // the traversal reaches the new layer
    this->SetVisited(N->GetID());
  };
  // traverse
  TraverseBreadth(GetRootNode(),repatch);
};


DAG::~DAG() {
  auto del = [] (auto m) {
    for(auto [k,v] : m)
      if(v) delete v;
    m.clear();
  };
  del(nodeMap);
  del(layerMap);
}
