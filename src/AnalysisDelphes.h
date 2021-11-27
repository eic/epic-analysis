#ifndef AnalysisDelphes_
#define AnalysisDelphes_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <stdexcept>

// root
#include "TChain.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TRegexp.h"

// delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

//#include "fastjet/contrib/Centauro.hh"
//#include "fastjet/plugins/Centauro/Centauro.hh"

// largex-eic
#include "Analysis.h"
#include "Histos.h"
#include "HistosDAG.h"
#include "Kinematics.h"
#include "CutDef.h"
#include "BinSet.h"
#include "SimpleTree.h"
#include "Weights.h"


class AnalysisDelphes : public Analysis
{
  public:
    AnalysisDelphes(
        TString infileName_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~AnalysisDelphes();

    // perform the analysis
    void Execute() override;

  private:

  ClassDefOverride(AnalysisDelphes,1);
};

#endif
