#include "DAGnode.h"

ClassImp(DAGnode)

using std::cout;
using std::cerr;
using std::endl;

// constructor
DAGnode::DAGnode(Int_t nodeType_, TString id_)
  : nodeType(nodeType_)
  , id(id_)
  , debug(true)
{
};

// print info
void DAGnode::Print() {
  TString nodeTypeStr;
  switch(nodeType) {
    case tBin: nodeTypeStr="Bin"; break;
    case tControl: nodeTypeStr="Control"; break;
    case tTop: nodeTypeStr="Top"; break;
    case tBottom: nodeTypeStr="Bottom"; break;
    default: nodeTypeStr="Unknown";
  };
  cout << "NODE ::: " << id << " ::: " << nodeTypeStr << endl;
  cout << "  Inputs:" << endl;
  for(DAGnode *N : inputList) cout << "    - " << N->GetID() << endl;
  cout << "  Outputs:" << endl;
  for(DAGnode *N : outputList) cout << "    - " << N->GetID() << endl;
};


// add node to given list
void DAGnode::AddInput(DAGnode *N) { AddNodeToList(N,inputList,"input"); };
void DAGnode::AddOutput(DAGnode *N) { AddNodeToList(N,outputList,"output"); };
void DAGnode::AddNodeToList(DAGnode *N, std::vector<DAGnode*> &list, TString listName) {
  if(debug) {
    cout << "ADD TO NODE: " << listName << " " << N->GetID() << " to " << this->GetID() << endl;
  };
  for(DAGnode *M : list) {
    if( M->GetID() == N->GetID() ) {
      cerr << "WARNING: tried to add duplicate node "
           << M->GetID() << "to " << listName << " list" << endl;
      return;
    };
  };
  list.push_back(N);
};


DAGnode::~DAGnode() {
};

