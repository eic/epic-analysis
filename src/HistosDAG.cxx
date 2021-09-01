#include "HistosDAG.h"

ClassImp(HistosDAG)

// default constructor
HistosDAG::HistosDAG()
  : debug(true)
{
  InitializeDAG();
};


// build the DAG from specified bin scheme
void HistosDAG::Build(std::map<TString,BinSet*> binSchemes_) {
  // copy binSchemes
  binSchemes = binSchemes_;
  // initialize DAG and histosMap
  InitializeDAG();
  histosMap.clear();
  // add one layer for each BinSet with nonzero bins defined
  for(auto kv : binSchemes) {
    if(debug) std::cout << "add BinSet " << kv.first << " to HistosDAG" << std::endl;
    BinSet *binScheme = kv.second;
    if(binScheme->GetNumBins()>0) AddLayer(binScheme);
  };
  // payload to create Histos objects
  Payload([this](NodePath P){
    TString histosN = "histos_";
    TString histosT = "";
    if(debug) std::cout << "At path " << Node::PathString(P) << ": ";
    // set name and title
    for(Node *N : P) {
      if(N->GetNodeType()==NT::bin) { // TODO: improve names and titles, sort them
        histosN += "_" + N->GetID();
        histosT += N->GetCut()->GetCutTitle() + ", ";
      };
    };
    if(debug) std::cout << "Create " << histosN << std::endl;
    // instantiate Histos object
    Histos *H = new Histos(histosN,histosT);
    // add CutDefs to Histos object
    for(Node *N : P) { if(N->GetNodeType()==NT::bin) H->AddCutDef(N->GetCut()); };
    // append to `histosMap`
    histosMap.insert(std::pair<NodePath,Histos*>(P,H));
  });
  // execution
  if(debug) std::cout << "Begin Histos instantiation..." << std::endl;
  ExecuteAndClearOps();
};


// payload wrapper operators, executed on the specified Histos object
// - lambda arguments: ( Histos* )
void HistosDAG::ForEach(std::function<void(Histos*)> op) {
  Payload( [op,this](NodePath P){ op(this->GetHistos(P)); } );
};
// - lambda arguments: ( Histos*, NodePath )
void HistosDAG::ForEach(std::function<void(Histos*,NodePath)> op) {
  Payload( [op,this](NodePath P){ op(this->GetHistos(P),P); } );
};


// return Histos* associated with the given NodePath
Histos *HistosDAG::GetHistos(NodePath P) {
  Histos *ret;
  try { ret = histosMap.at(P); }
  catch(const std::out_of_range &ex) {
    std::cerr << "ERROR: no Histos associated with NodePath "
              << Node::PathString(P) << std::endl;
    return nullptr;
  };
  return ret;
};


HistosDAG::~HistosDAG() {
};

