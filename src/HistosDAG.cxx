#include "HistosDAG.h"

ClassImp(HistosDAG)

// default constructor
HistosDAG::HistosDAG()
  : debug(false)
{
  InitializeDAG();
};


// build the DAG from specified bin scheme
void HistosDAG::Build(std::map<TString,BinSet*> binSchemes) {
  // initialize DAG and histosMap
  InitializeDAG();
  histosMap.clear();
  // add the finalState layer first, if it exists
  try { 
    BinSet *finalLayer = binSchemes.at("finalState");
    if(finalLayer->GetNumBins()>0) AddLayer(finalLayer);
  } catch(const std::out_of_range &ex) {
    std::cerr << "WARNING: no finalState bins defined" << std::endl;
  };
  // add one layer for each BinSet with nonzero bins defined
  for(auto kv : binSchemes) {
    if(debug) std::cout << "add BinSet " << kv.first << " to HistosDAG" << std::endl;
    if(kv.first=="finalState") continue;
    BinSet *binScheme = kv.second;
    if(binScheme->GetNumBins()>0) AddLayer(binScheme);
  };
  // leaf operator, to create Histos objects
  LeafOp([this](NodePath *P){
    TString histosN = "histos";
    TString histosT = "";
    if(debug) std::cout << "At path " << P->PathString() << ": ";
    // set name and title
    for(Node *N : P->GetSortedBins()) {
      histosN += "__" + N->GetID();
      histosT += N->GetCut()->GetCutTitle() + ", ";
    };
    if(debug) std::cout << "Create " << histosN << std::endl;
    // instantiate Histos object
    Histos *H = new Histos(histosN,histosT);
    // add CutDefs to Histos object
    for(Node *N : P->GetBinNodes()) { H->AddCutDef(N->GetCut()); };
    // append to `histosMap`
    histosMap.insert(std::pair<std::set<Node*>,Histos*>(P->GetBinNodes(),H));
  });
  // execution
  if(debug) std::cout << "Begin Histos instantiation..." << std::endl;
  ExecuteAndClearOps();
};


// build the DAG from ROOT file; all BinSets will become layers and
// all Histos objects will be linked to NodePaths
void HistosDAG::Build(TFile *rootFile) {
  // initialize DAG and histosMap, and read rootFile keys
  InitializeDAG();
  histosMap.clear();
  TListIter nextKey(rootFile->GetListOfKeys());
  TString keyname;
  // add each BinSet as a new layer
  while(TKey *key = (TKey*)nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp("^binset__"))) {
      if(debug) std::cout << "READ LAYER " << keyname << std::endl;
      BinSet *B = (BinSet*)key->ReadObj();
      if(B->GetNumBins()>0) AddLayer(B);
    };
  };
  nextKey.Reset();
  // add each Histos to histMap
  while(TKey *key = (TKey*)nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp("^histos__"))) {
      // get NodePath from Histos name
      if(debug) std::cout << "READ HISTOS " << keyname << std::endl;
      NodePath P;
      P.nodes.insert(GetRootNode());
      P.nodes.insert(GetLeafNode());
      TString tokID;
      Ssiz_t tf=0;
      while(keyname.Tokenize(tokID,tf,"__")) {
        if(tokID=="histos") continue;
        Node *N = GetNode(tokID);
        if(N) P.nodes.insert(N);
        else {
          std::cerr << "ERROR: mismatch of Node \"" << tokID << "\" between Histos and BinSets" << std::endl;
          return;
        };
      };
      // append to `histosMap`
      if(debug) std::cout << "-> PATH: " << P.PathString() << std::endl;
      histosMap.insert(std::pair<std::set<Node*>,Histos*>(P.GetBinNodes(),(Histos*)key->ReadObj()));
    };
  };
};


// return Histos* associated with the given NodePath
Histos *HistosDAG::GetHistos(NodePath *P) {
  Histos *ret;
  try { ret = histosMap.at(P->GetBinNodes()); }
  catch(const std::out_of_range &ex) {
    std::cerr << "ERROR: no Histos associated with NodePath "
              << P->PathString() << std::endl;
    return nullptr;
  };
  return ret;
};

// return Histos* associated with the given external NodePath, by ID-matching its Nodes to the local DAG's Nodes
Histos *HistosDAG::GetHistosExternal(NodePath *extP) {
  NodePath *intP = new NodePath();
  for(auto extN : extP->nodes) {
    auto intN = this->GetNode(extN->GetID(),true); // find internal node by ID-matching external node
    if(intN!=nullptr) intP->nodes.insert(intN);
  }
  return this->GetHistos(intP);
};

HistosDAG::~HistosDAG() {
};

