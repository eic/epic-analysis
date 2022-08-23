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
