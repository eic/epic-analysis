/* ADAGE front-end template class
 * A - Analysis in a
 * D - Directed
 * A - Acyclic
 * G - Graph
 * E - Environment
 *
 * - Associates DAG paths to unique objects
 * - User-friendly higher order functions (`Payload`, `MultiPayload`, ...) for operating
 *   on the stored objects
 * - Inherits from back-end class `DAG`
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include <TString.h>
#include <TRegexp.h>
#include <TFile.h>
#include <TKey.h>

// adage
#include "DAG.h"
#include "BinSet.h"
#include "Node.h"
#include "NodePath.h"


template<class PL> class Adage : public DAG
{
  public:

    Adage(TString payloadNamePrefix_="") : DAG() {
      payloadNamePrefix = payloadNamePrefix_; // optional unique name prefix, used only for storage in TFiles
    };
    ~Adage() {};

    
    // METHODS
    // -----------------------------------------------------------------------------------

    // virtual method to build the DAG and instantiate payload data, from specified bin scheme
    /* - derived classes must override this, following this general format:
     *   - first, call `BuildDAG(binSchemes)`
     *   - next, stage a `LeafOp` lambda:
     *     - instantiate your payload object, give it a unique name (`CreatePayloadName`, `CreatePayloadTitle`), etc.
     *     - the last line of the `LeafOp` lambda should be `InsertPayloadData`
     *   - last, call `ExecuteAndClearOps()`;
     */
    virtual void Build(std::map<TString,BinSet*> binSchemes) = 0;

    // build the DAG from a ROOT file
    /* - all BinSets will become layers
     * - all objects that have names starting with `payloadNamePrefix` will be associated to DAG paths
     */
    void BuildFromFile(TFile *rootFile);

    // return payload data (`PL` pointer) associated with the given `NodePath`
    /* - if you have a `NodePath` from another DAG that has the same binning
     *    scheme, set `NodePath_is_external=true`
     */
    PL *GetPayloadData(NodePath *P, bool NodePath_is_external=false);


    // OPERATORS
    // -----------------------------------------------------------------------------------

    // payload operator, executed on the specified object; see `FormatPayload`
    // for allowed arguments of operator `op`
    template<class O>
      void Payload(O op) {
        auto opFormatted = FormatPayload(op);
        LeafOp( [opFormatted,this](NodePath *P){ opFormatted(this->GetPayloadData(P),P); } );
      };

    // multi-payload operators, for defining multiple payloads for a subloop
    template<class O1, class O2, class O3>
      void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore, O3 opAfter) {
        auto opFormatted = FormatPayload(opPayload);
        MultiLeafOp(
            layers,
            [opFormatted,this](NodePath *P){ opFormatted(this->GetPayloadData(P),P); },
            opBefore,
            opAfter
            );
      };
    template<class O1, class O2>
      void MultiPayload(std::vector<TString> layers, O1 opPayload, O2 opBefore) { MultiPayload(layers,opPayload,opBefore,[](){}); };
    template<class O1>
      void MultiPayload(std::vector<TString> layers, O1 opPayload) { MultiPayload(layers,opPayload,[](){},[](){}); };

     
  protected: // -----------------------------------------------------------------------------------

    // build the DAG from specified bin scheme
    void BuildDAG(std::map<TString,BinSet*> binSchemes, std::vector<TString> firstLayers={});

    // create a unique name or title for a payload object
    TString CreatePayloadName(NodePath *P);
    TString CreatePayloadTitle(TString titlePrefix, NodePath *P);

    // associate a path to payload data
    void InsertPayloadData(NodePath *P, PL *PLptr);


  private: // -----------------------------------------------------------------------------------

    // format payload operators with proper arguments `(PayLoad object pointer, NodePath pointer)`, to allow easy overloading
    static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(PL*,NodePath*)> op) { return op;     };
    static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(NodePath*,PL*)> op) { return [op](PL *H, NodePath *P){ op(P,H); }; };
    static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(PL*)>           op) { return [op](PL *H, NodePath *P){ op(H);   }; };
    static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void(NodePath*)>     op) { return [op](PL *H, NodePath *P){ op(P);   }; };
    static std::function<void(PL*,NodePath*)> FormatPayload(std::function<void()>              op) { return [op](PL *H, NodePath *P){ op();    }; };

    TString payloadNamePrefix;                 // unique name, used for TFile I/O
    std::map<std::set<Node*>,PL*> payloadHash; // map DAG path of bin nodes -> PL*

  ClassDefOverride(Adage,1);
};


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


// build the DAG from specified bin scheme
template <class PL>
void Adage<PL>::BuildDAG(std::map<TString,BinSet*> binSchemes, std::vector<TString> firstLayers) {
  // initialize DAG and payloadHash
  InitializeDAG();
  payloadHash.clear();
  // add `firstLayers` first, if any, and if they exist in `binSchemes`
  for(TString firstLayerN : firstLayers) {
    try { 
      if(debug) std::cout << "add BinSet " << firstLayerN << " to Adage" << std::endl;
      BinSet *firstLayer = binSchemes.at(firstLayerN);
      if(firstLayer->GetNumBins()>0) AddLayer(firstLayer);
    } catch(const std::out_of_range &ex) {
      std::cerr << "WARNING: no " << firstLayerN << " bins defined" << std::endl;
    };
  };
  // add one layer for each BinSet with nonzero bins defined
  for(auto kv : binSchemes) {
    if( std::find(firstLayers.begin(), firstLayers.end(), kv.first) != firstLayers.end() ) continue;
    if(debug) std::cout << "add BinSet " << kv.first << " to Adage" << std::endl;
    BinSet *binScheme = kv.second;
    if(binScheme->GetNumBins()>0) AddLayer(binScheme);
  };
};

// build the DAG from ROOT file
template <class PL>
void Adage<PL>::BuildFromFile(TFile *rootFile) {
  // initialize DAG and payloadHash, and read rootFile keys
  InitializeDAG();
  payloadHash.clear();
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
  // add each payload found to `payloadHash`
  while(TKey *key = (TKey*)nextKey()) {
    keyname = TString(key->GetName());
    if(keyname.Contains(TRegexp(TString("^")+payloadNamePrefix))) {
      // get NodePath from matched name
      if(debug) std::cout << "READ " << payloadNamePrefix << " " << keyname << std::endl;
      NodePath P;
      P.nodes.insert(GetRootNode());
      P.nodes.insert(GetLeafNode());
      TString tokID;
      Ssiz_t tf=0;
      while(keyname.Tokenize(tokID,tf,"__")) {
        if(tokID==payloadNamePrefix) continue;
        Node *N = GetNode(tokID);
        if(N) P.nodes.insert(N);
        else {
          std::cerr << "ERROR: mismatch of Node \"" << tokID << "\" between Payload Object and BinSets" << std::endl;
          return;
        };
      };
      // append to `payloadHash`
      if(debug) std::cout << "-> PATH: " << P.PathString() << std::endl;
      payloadHash.insert({P.GetBinNodes(),(PL*)key->ReadObj()});
    };
  };
};

// create a unique name for the payload object associated to NodePath `P`; it will begin with `payloadNamePrefix`
template <class PL>
TString Adage<PL>::CreatePayloadName(NodePath *P) {
  TString ret = payloadNamePrefix;
  for(Node *N : P->GetSortedBins()) {
    ret += "__" + N->GetID();
  };
  return ret;
};

// create a title for the payload object associated to NodePath `P`; it will include CutDef titles
template <class PL>
TString Adage<PL>::CreatePayloadTitle(TString titlePrefix, NodePath *P) {
  TString ret = titlePrefix;
  for(Node *N : P->GetSortedBins()) {
    ret += N->GetCut()->GetCutTitle() + ", ";
  };
  ret(TRegexp(", $")) = "";
  return ret;
};

// add payload object to `payloadHash`, which associates a NodePath to that object
template <class PL>
void Adage<PL>::InsertPayloadData(NodePath *P, PL *PLptr) {
  payloadHash.insert({P->GetBinNodes(),PLptr});
};

// return payload object associated with the given NodePath
template <class PL>
PL *Adage<PL>::GetPayloadData(NodePath *P,  bool NodePath_is_external) {
  PL *ret;
  NodePath *intP;
  if(NodePath_is_external) {
    intP = new NodePath();
    for(auto extN : P->nodes) {
      auto intN = this->GetNode(extN->GetID(),true); // find internal node by ID-matching external node
      if(intN!=nullptr) intP->nodes.insert(intN);
    }
  } else intP=P;
  try { ret = payloadHash.at(intP->GetBinNodes()); }
  catch(const std::out_of_range &ex) {
    std::cerr << "ERROR: no Payload Object associated with NodePath "
              << intP->PathString() << std::endl;
    return nullptr;
  };
  return ret;
};
