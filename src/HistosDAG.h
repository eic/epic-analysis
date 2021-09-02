#ifndef HistosDAG_
#define HistosDAG_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TKey.h"

// largex-eic
#include "DAG.h"
#include "Histos.h"
#include "BinSet.h"


class HistosDAG : public DAG
{
  public:
    HistosDAG();
    ~HistosDAG();

    // build the DAG from specified bin scheme
    void Build(std::map<TString,BinSet*> binSchemes);

    // build the DAG from ROOT file; all BinSets will become layers and
    // all Histos objects will be linked to NodePaths
    void Build(TFile *rootFile);

    // payload operator, executed on the specified Histos object
    void ForEach(std::function<void(Histos*)> op);
    void ForEach(std::function<void(Histos*,NodePath)> op);

    // return Histos* associated with the given NodePath
    Histos *GetHistos(NodePath P);

  private:
    Bool_t debug;
    std::map<NodePath,Histos*> histosMap; // map full DAG path -> Histos*

  ClassDefOverride(HistosDAG,1);
};

#endif
