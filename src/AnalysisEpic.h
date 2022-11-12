#ifndef AnalysisEpic_
#define AnalysisEpic_

// data model
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "podio/CollectionBase.h"
#include "edm4hep/utils/kinematics.h"

// data model collections
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/InclusiveKinematicsCollection.h"

// sidis-eic
#include "Analysis.h"


class AnalysisEpic : public Analysis
{
  public:
    AnalysisEpic(TString infileName_="", TString outfilePrefix_="");
    ~AnalysisEpic();

    void Execute() override;

    // settings
    bool crossCheckKinematics;

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
    void PrintParticle(const edm4eic::ReconstructedParticle& P);
    void PrintAssociatedParticles(
        const edm4hep::MCParticle& simPart,
        const edm4eic::ReconstructedParticle& recPart
        );


  protected:

    // get PDG from reconstructed particle
    int GetReconstructedPDG(
        const edm4hep::MCParticle& simPart,
        const edm4eic::ReconstructedParticle& recPart
        );
    // run `payload` for all [Reconstructed Particle] <-> [MC Particle] associations
    // payload signature: (simPart, recPart, reconstructed PDG)
    void LoopMCRecoAssocs(
        const edm4eic::MCRecoParticleAssociationCollection& mcRecAssocs,
        std::function<void(const edm4hep::MCParticle&, const edm4eic::ReconstructedParticle&, int)> payload,
        bool printParticles=false
        );

  private:
    podio::ROOTReader podioReader;
    podio::EventStore evStore;

    // reconstructed PDG cache table
    bool useCachedPDG;
    std::map<std::pair<unsigned int,unsigned int>, int> pdgCache; // map : {simPart.id(),recPart.id()} -> recPDG

    ClassDefOverride(AnalysisEpic,1);
};

#endif
