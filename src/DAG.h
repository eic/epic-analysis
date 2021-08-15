#ifndef DAG_
#define DAG_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <functional>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"


// largex-eic
#include "CutDef.h"
#include "Node.h"

class DAG : public TObject
{
  public:

    DAG();
    ~DAG();

    // initialize DAG, with only the root node connected to the leaf node
    void Initialize();

    // node accessors
    // - search for a node by ID; return nullptr if not found;
    //   set silence=true for no error print out when node not found
    Node *GetNode(TString id_, Bool_t silence=false);
    // get root, leaf, or any other unique node; non-uniqueness 
    // will cause error printout
    Node *GetRootNode();
    Node *GetLeafNode();
    Node *GetUniqueNode(Int_t type_,TString typeStr);

    // add nodes and edges, which are directional connectors between nodes
    // - nodes are only added if the ID is unique; if not unique,
    //   do nothing (silence=true suppresses error print out)
    void AddNode(Node *N, Bool_t silence=false);
    void AddNode(Int_t nodeType, TString id, Bool_t silence=false);
    void AddEdge(Node *inN, Node *outN);
    void ModifyNode(Node *N, TString newName, Int_t newType=-1); // rename or repurpose a node

    // add a layer of nodes, fully connected to the last layer of the DAG
    // - primary usage is to add a layer of bins from a BinSet
    void AddLayer(std::vector<Node*> nodes);

    // DAG traversals: iterate through nodes, executing the lambda on each;
    // breadth-first and depth-first traversals are supported
    void TraverseBreadth(Node *N, std::function<void(Node*)> lambda);
    void TraverseDepth(Node *N, std::function<void(Node*)> lambda);

    // simplify DAG: convert all control nodes into full connections between
    // the adjacent layers; all control nodes will be removed; you can also
    // convert a single control node to a full connection
    void Simplify();
    void RemoveControl(TString id);
    void RemoveControl(Node *N);

    // remove nodes and edges
    void RemoveNode(Node *N);
    void RemoveEdge(Node *inN, Node *outN);

    // print the whole DAG (breadth first)
    void Print(TString header="DAG");


  protected:
    Bool_t Visited(TString id_);

  private:
    Bool_t debug;
    std::map<TString,Node*> nodeMap;
    std::vector<TString> visitList;

  ClassDef(DAG,1);
};

#endif
