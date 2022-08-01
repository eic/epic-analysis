// Implement Histos payload for Adage
#ifndef HistosDAG_
#define HistosDAG_

#include "Adage.h"
#include "Histos.h"

class HistosDAG : public Adage<Histos>
{
  public:
    HistosDAG();
    ~HistosDAG();

    // build the DAG from specified bin scheme
    void Build(std::map<TString,BinSet*> binSchemes) override;

  ClassDefOverride(HistosDAG,1);
};

#endif
