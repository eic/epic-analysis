/* ADAGE backend data structure: a Directed Acyclic Graph (DAG)
 * A - Analysis in a
 * D - Directed
 * A - Acyclic
 * G - Graph
 * E - Environment
 */
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
#include "TObjArray.h"
#include "TNamed.h"
#include "TString.h"


// largex-eic
#include "CutDef.h"
#include "BinSet.h"
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
    void AddEdge(Node *inN, Node *outN, Bool_t silence=false);
    void ModifyNode(Node *N, TString newName, Int_t newType=-1); // rename or repurpose a node

    // add a layer of nodes, fully connected ("patched") to the last layer of the DAG
    // - primary usage is to add a layer of bins from a BinSet
    // - "patch" refers to any set of edges between two layers
    void AddLayer(BinSet *BS);
    void AddLayer(std::vector<Node*> nodes);

    // DAG traversals: iterate through nodes, executing the lambda on each;
    // breadth-first and depth-first traversals are supported
    void TraverseBreadth(Node *N, std::function<void(Node*)> lambda);
    void TraverseDepth(Node *N, std::function<void(Node*)> lambda);

    // patch operations: manipulate connections ("patches") between layers
    /*   - "patch" refers to a set of edges between two layers, which
     *     can be two types:
     *     - "full patch": each node of one layer is connected to each node of the other
     *     - "control patch": each node of one layer is connected to a control node, which
     *       is then connected to each node of the next layer
     */
    // convert control patch(es) to full patch(es); control node(s) will be removed
    void RepatchToFull(TString id);
    void RepatchToFull(Node *N);
    void RepatchAllToFull();
    // re-patch a layer to the current leaf node's output, and re-patch the
    // adjacent layers together; the current leaf node will become a control
    // node, and a new leaf node is created
    void RepatchToLeaf(TString varName);
    // mark a node as visited
    void SetVisited(TString id_) { visitList.push_back(id_); };

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
