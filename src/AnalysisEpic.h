#ifndef AnalysisEpic_
#define AnalysisEpic_

// data model
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/utils/kinematics.h"
#include "eicd/InclusiveKinematicsCollection.h"

// sidis-eic
#include "Analysis.h"


class AnalysisEpic : public Analysis
{
  public:
    AnalysisEpic(
        TString infileName_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~AnalysisEpic();

    void Execute() override;

  private:
    podio::ROOTReader podioReader;
    podio::EventStore evStore;

    ClassDefOverride(AnalysisEpic,1);
};

#endif
