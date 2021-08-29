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

// largex-eic
#include "DAG.h"
#include "Histos.h"


class HistosDAG : public DAG
{
  public:
    HistosDAG();
    ~HistosDAG();

    /* operator staging templates:
     * - `op` should be a lambda, which is allowed to have
     *   arguments `Node*` or `NodePath`, or both, or neither
     */

    // initial operator, executed at the beginning of DAG::ExecuteOps
    template<class O> void Initial(O op) { this->GetRootNode()->StageInboundOp(op); };

    // final operator, executed at the end of DAG::ExecuteOps
    template<class O> void Final(O op) { this->GetRootNode()->StageOutboundOp(op); };


    // combination of Before and After
    template<class O1, class O2> void Control(/*TString layer*/ std::vector<TString> layers, O1 opBefore, O2 opAfter) {
      if(layers.size()==0) {
        std::cerr << "ERROR: empty layers list in HistosDAG::Control" << std::endl;
        return;
      };
      bool first = true;
      TString controlID;
      for(TString layer : layers) {
        RepatchToLeaf(layer);
        if(first) controlID = layer+"__control";
        else RepatchToFull(layer+"__control");
        first = false;
      };
      GetNode(controlID)->StageInboundOp(opBefore);
      GetNode(controlID)->StageOutboundOp(opBefore);
    };


    // control operator, executed before the listed layers; these layers will be moved to the end
    template<class O> void Before(std::vector<TString> layers, O op) {
      // TODO
      // also be careful Before, After, and Control don't all overwrite each other
    };


    // control operator, executed after the listed layers; these layers will be moved to the end
    template<class O> void After(std::vector<TString> layers, O op) {
      // TODO
      // also be careful Before, After, and Control don't all overwrite each other
    };


    // payload operator, executed on the specified Histos object
    void Payload(std::function<void(Histos*)> op) {
      // TODO: implement this; currently just test code
      GetLeafNode()->StageInboundOp([](NodePath P){
        std::cout << "ALGORITHM on ";
        Node::PrintPath(P);
      });
    };


    // traverse DAG and run all operators
    void Execute() { ExecuteOps(); };


  private:
    Bool_t debug;

  ClassDefOverride(HistosDAG,1);
};

#endif
