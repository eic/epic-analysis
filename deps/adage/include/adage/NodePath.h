#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include <vector>

// ROOT
#include <TSystem.h>
#include <TObject.h>
#include <TNamed.h>
#include <TString.h>
#include <TRegexp.h>

// adage
#include "CutDef.h"
#include "Node.h"

class Node;

class NodePath : public TObject
{
  public:

    NodePath();
    ~NodePath();

    // user-friendly methods
    /* - NodePath is effectively a linked list
     * - use these methods to get specific nodes, along with things the nodes hold such as cuts
     * - this class stores the list unsorted; use the sort methods to traverse the linked list in order 
     */
    TString BinListName(); // return a string of bin names, formatted for object names: `"bin1_bin2_bin3"`
    TString BinListString(); // return a string of bin names, formatted for printing: `"[ bin1 bin2 bin3 ]"`
    void PrintBinList(); // print a string of node names
    TString CutListString(); // return a string of cuts (CutDef titles)
    Node *GetBinNode(TString varName); // return the bin node for the given variable name
    Node *GetPreviousNode(Node *N, bool silence=false); // return the node that has `N` as one of its outputs
    Node *GetNextNode(Node *N, bool silence=false); // return the node that has `N` as one of its inputs
    Node *GetFirstNode(); // return the first node in the path
    Node *GetLastNode(); // return the first node in the path
    std::set<Node*> GetBinNodes(); // get the full set of bin nodes, useful for iteration; if you need all 
                                   // nodes (including leaf, root, and control nodes), use `NodePath::nodes`)
    std::vector<Node*> GetSortedPath(); // returns path in DAG order (NodePath encapsulates an unordered set of nodes)
    std::vector<Node*> GetSortedBins(); // returns list of bins in DAG order

    // additional methods
    TString PathString(); // return a string of node names
    void PrintPath(); // print a string of node names

    // encapsulation of std::set<Node*>
    // - this is the internal data structure of NodePath; it is public to grant
    //   full access to std::set functionality
    // - it is an unordered set, since we want two NodePaths which contain the same set
    //   of nodes to be equivalent, even if the DAG connections are not (cf. `GetSortedPath()`)
    std::set<Node*> nodes;

  private:

  ClassDef(NodePath,1);
};
