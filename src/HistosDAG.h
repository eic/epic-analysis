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
#include "Node.h"
#include "NodePath.h"


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

    // payload operator, executed on the specified Histos object; see `FormatPayload`
    // for allowed arguments of operator `op`
    template<class O>
    void Payload(O op) {
      auto opFormatted = FormatPayload(op);
      LeafOp( [opFormatted,this](NodePath *P){ opFormatted(this->GetHistos(P),P); } );
    };

    // multi-payload operators, for defining multiple payloads for a subloop
    template<class O1, class O2, class O3>
    void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore, O3 opAfter) {
      auto opFormatted = FormatPayload(opPayload);
      MultiLeafOp(
          layers,
          [opFormatted,this](NodePath *P){ opFormatted(this->GetHistos(P),P); },
          opBefore,
          opAfter
          );
    };
    template<class O1, class O2>
    void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore) { MultiPayload(layers,opPayload,opBefore,[](){}); };
    template<class O1>
    void MultiPayload(std::vector<TString> layers, O1 opPayload) { MultiPayload(layers,opPayload,[](){},[](){}); };

    // format payload operators with proper arguments, to allow easy overloading
    static std::function<void(Histos*,NodePath*)> FormatPayload(std::function<void(Histos*,NodePath*)> op) {
      return op;
    };
    static std::function<void(Histos*,NodePath*)> FormatPayload(std::function<void(NodePath*,Histos*)> op) {
      return [op](Histos *H, NodePath *P){ op(P,H); };
    };
    static std::function<void(Histos*,NodePath*)> FormatPayload(std::function<void(Histos*)> op) {
      return [op](Histos *H, NodePath *P){ op(H); };
    };

    // return Histos* associated with the given NodePath
    Histos *GetHistos(NodePath *P);

  private:
    Bool_t debug;
    std::map<std::set<Node*>,Histos*> histosMap; // map DAG path of bin nodes -> Histos*

  ClassDefOverride(HistosDAG,1);
};

#endif
