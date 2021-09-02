#ifndef NodePath_
#define NodePath_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"

// largex-eic
#include "CutDef.h"
#include "Node.h"

class Node;

class NodePath : public TObject
{
  public:

    NodePath();
    ~NodePath();

    // user-friendly methods
    TString BinListString(); // return a string of bin names
    void PrintBinList(); // print a string of node names
    TString CutListString(); // return a string of cuts (CutDef titles)
    // get the full set of bin nodes, useful for iteration; if you need
    // all nodes (including leaf, root, and control nodes), use `NodePath::nodes`)
    std::set<Node*> GetBinNodes();

    // additional methods
    TString PathString(); // return a string of node names
    void PrintPath(); // print a string of node names

    // encapsulation of std::set<Node*>
    // - this is the internal data structure of NodePath; it is public to grant
    //   full access to std::set functionality
    std::set<Node*> nodes;

  private:

  ClassDef(NodePath,1);
};

#endif
