// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Kevin Adkins, Brian Page, Connor Pecar
#ifndef EXCLUDE_DELPHES

/* NOTE:
 * if you make changes, MAINTAIN DOCUMENTATION IN ../doc/kinematicsJets.md
 */

#pragma once

// Kinematics base calss
#include "Kinematics.h"

// Delphes includes
#include <classes/DelphesClasses.h>
#include <fastjet/ClusterSequence.hh>
#ifdef INCLUDE_CENTAURO
#include <fastjet/plugins/Centauro/Centauro.hh>
#endif // CENTAURO

using std::map;
using std::cout;
using std::cerr;
using std::endl;

class KinematicsJets : public Kinematics
{
  public:
  KinematicsJets(Double_t enEleBeam, Double_t enIonBeam, Double_t crossAng);
    ~KinematicsJets();

    // DELPHES-specific methods ////////////////////////// -- to be removed when jets are added to AnalysisEpic

    // jet calculators
    void GetJets(
        TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle,
	int jetAlgo, double jetRadius, double jetMinPt
        );
    void CalculateJetKinematics(fastjet::PseudoJet jet);
    void CalculateJetResolution(double deltaRCut);
#ifdef INCLUDE_CENTAURO
    void GetBreitFrameJets(
        TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
        );
    void CalculateBreitJetKinematics(fastjet::PseudoJet jet);
#endif // ifdef INCLUDE_CENTAURO

    // end DELPHES-specific methods //////////////////////////
  
    // jet objects
    int jetAlgo;
    double jetRad, jetMinPt;

    std::vector<fastjet::PseudoJet> jetsRec, jetsTrue;
    std::vector<fastjet::PseudoJet> breitJetsRec, breitJetsTrue;
    std::map<double, int> jetConstituents;
    fastjet::ClusterSequence csRec;
    fastjet::ClusterSequence csTrue;

    Double_t zjet, pTjet, qTjet, mTjet, etajet, phijet, mjet, ejet;
    Double_t deltaRjet;
    int matchStatusjet;
    Double_t pTmtjet, mTmtjet, etamtjet, phimtjet, mmtjet, emtjet;
    std::vector<double> jperp;
    std::vector<double> zhad_jet;
    // struck quark information
    Double_t quarkpT;

  ClassDef(KinematicsJets,1);
};

#endif
