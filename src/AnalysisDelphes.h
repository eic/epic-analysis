// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks, Connor Pecar, Matthew McEneaney

#ifndef AnalysisDelphes_
#define AnalysisDelphes_

// delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

//#include "fastjet/contrib/Centauro.hh"
//#include "fastjet/plugins/Centauro/Centauro.hh"

// sidis-eic
#include "Analysis.h"


class AnalysisDelphes : public Analysis
{
  public:
    AnalysisDelphes(TString configFileName_="", TString outfilePrefix_="");
    ~AnalysisDelphes();

    // perform the analysis
    void Execute() override;

  private:

  ClassDefOverride(AnalysisDelphes,1);
};

#endif
