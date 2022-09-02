#ifndef AnalysisEpic_
#define AnalysisEpic_

// data model
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "podio/CollectionBase.h"
#include "edm4hep/utils/kinematics.h"

// data model collections
#include "edm4hep/MCParticleCollection.h"
#include "eicd/ReconstructedParticleCollection.h"
#include "eicd/MCRecoParticleAssociationCollection.h"
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

    // return Lorentz vector for a given particle
    template <class ParticleType>
    TLorentzVector GetP4(ParticleType& P) {
      return TLorentzVector(
          P.getMomentum().x,
          P.getMomentum().y,
          P.getMomentum().z,
          P.getEnergy()
          );
    }

    // printers
    void PrintParticle(const edm4hep::MCParticle& P);
    void PrintParticle(const eicd::ReconstructedParticle& P);


  protected:

    // get PDG from reconstructed particle
    int GetReconstructedPDG(
        const edm4hep::MCParticle& simPart,
        const eicd::ReconstructedParticle& recPart,
        bool& usedTruth
        );

  private:
    podio::ROOTReader podioReader;
    podio::EventStore evStore;

    ClassDefOverride(AnalysisEpic,1);
};

#endif
