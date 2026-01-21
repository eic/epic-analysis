// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

#pragma once

// PODIO
#include <podio/ROOTFrameReader.h>
#include <podio/Frame.h>

// data model
#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
#include <edm4eic/MCRecoParticleAssociationCollection.h>
#include <edm4eic/InclusiveKinematicsCollection.h>

// utilities
#include <edm4hep/utils/kinematics.h>

// epic-analysis
#include "Analysis.h"


class AnalysisEpicPodio : public Analysis
{
  public:
    AnalysisEpicPodio(TString infileName_="", TString outfilePrefix_="");
    ~AnalysisEpicPodio();

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


  protected:

    // get PDG from reconstructed particle; resort to true PDG, if
    // PID is unavailable (sets `usedTruth` to true)
    int GetPDG(
        const edm4hep::MCParticle& simPart,
        const edm4eic::ReconstructedParticle& recPart,
        bool& usedTruth
        );
    // common loop over Reconstructed Particle <-> MC Particle associations
    // payload signature: (simPart, recPart, reconstructed PDG)
    void LoopMCRecoAssocs(
        const edm4eic::ReconstructedParticleCollection&     recParts,
        const edm4eic::MCRecoParticleAssociationCollection& mcRecAssocs,
        std::function<void(const edm4hep::MCParticle&, const edm4eic::ReconstructedParticle&, int)> payload,
        bool printParticles=false
        );

  private:
    podio::ROOTFrameReader podioReader;

    ClassDefOverride(AnalysisEpicPodio,1);
};
