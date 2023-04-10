// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Kevin Adkins, Brian Page, Connor Pecar
#ifndef EXCLUDE_DELPHES

/* NOTE:
 * if you make changes, MAINTAIN DOCUMENTATION IN ../doc/kinematicsJets.md
 */

#include "KinematicsJets.h"
ClassImp(KinematicsJets)

KinematicsJets::KinematicsJets(Double_t enEleBeam_, Double_t enIonBeam_, Double_t crossAng_) : Kinematics(enEleBeam_,enIonBeam_,crossAng_)
{/*Note: Energy units are GeV, angle units are mrad*/};


// JET METHODS ///////////////////////

void KinematicsJets::GetJets(
    TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
    TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle,
    int jetAlgo, double jetRadius, double jetMinPt
    )
{
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();
  itParticle.Reset();

  while(GenParticle *partTrue = (GenParticle*)itParticle() ){
    if( (partTrue->PID == 1 || partTrue->PID == 2) && (partTrue->Status == 23) ){
      // Status: 23->outgoing, but there's also 63->outgoing beam remnant. TODO: Which do we want?
      // from pythia 8 documentation
      quarkpT = partTrue->PT;
    }
  }

  std::vector<fastjet::PseudoJet> particles;
  std::vector<fastjet::PseudoJet> particlesTrue;
  jetConstituents.clear();
  // looping over final state particles, adding to particles vector
  while(Track *eflowTrack = (Track*)itEFlowTrack() ){
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
    int partPID = std::abs(trackParticle->PID);
    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) < 4.0 && eflowTrack->PT > 0.1){
        this->TransformToHeadOnFrame(eflowTrackp4,eflowTrackp4);
        if(partPID != 11) particles.push_back(fastjet::PseudoJet(eflowTrackp4.Px(),eflowTrackp4.Py(),eflowTrackp4.Pz(),eflowTrackp4.E()));

        //GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
        TLorentzVector partp4 = trackParticle->P4();
        this->TransformToHeadOnFrame(partp4,partp4);
        if(partPID != 11) particlesTrue.push_back(fastjet::PseudoJet(partp4.Px(),partp4.Py(),partp4.Pz(),partp4.E()));

        jetConstituents.insert({eflowTrackp4.Px(), eflowTrack->PID});
      }
    }
  }
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0){
        this->TransformToHeadOnFrame(towerPhotonp4,towerPhotonp4);
        particles.push_back(fastjet::PseudoJet(towerPhotonp4.Px(),towerPhotonp4.Py(),towerPhotonp4.Pz(),towerPhotonp4.E()));

        for(int i = 0; i < towerPhoton->Particles.GetEntries(); i++){
          GenParticle *photonPart = (GenParticle*)towerPhoton->Particles.At(i);
          TLorentzVector photonp4 = photonPart->P4();
          this->TransformToHeadOnFrame(photonp4,photonp4);
          particlesTrue.push_back(fastjet::PseudoJet(photonp4.Px(),photonp4.Py(),photonp4.Pz(),photonp4.E()));
        }
      }
    }
  }

  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if(!isnan(towerNeutralHadronp4.E())){
      if( std::abs(towerNeutralHadron->Eta) < 4.0){
        this->TransformToHeadOnFrame(towerNeutralHadronp4,towerNeutralHadronp4);
        particles.push_back(
          fastjet::PseudoJet(towerNeutralHadronp4.Px(),towerNeutralHadronp4.Py(),towerNeutralHadronp4.Pz(),towerNeutralHadronp4.E())
          );

        for(int i = 0; i < towerNeutralHadron->Particles.GetEntries(); i++){
          GenParticle *nhadPart = (GenParticle*)towerNeutralHadron->Particles.At(i);
          TLorentzVector nhadp4 = nhadPart->P4();
          this->TransformToHeadOnFrame(nhadp4,nhadp4);
          particlesTrue.push_back(fastjet::PseudoJet(nhadp4.Px(),nhadp4.Py(),nhadp4.Pz(),nhadp4.E()));
        }
      }
    }
  }

  // Set Jet Definition
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, jetRadius);

  // Redefine Jet Algo Based on Macro Input
  // Default is jetAlgo = 0 which is kt_algorithm
  if(jetAlgo == 1) jet_def.set_jet_algorithm(fastjet::cambridge_algorithm);
  if(jetAlgo == 2) jet_def.set_jet_algorithm(fastjet::antikt_algorithm);

  csRec = fastjet::ClusterSequence(particles, jet_def);
  csTrue = fastjet::ClusterSequence(particlesTrue, jet_def);
  jetsRec = sorted_by_pt(csRec.inclusive_jets(jetMinPt));
  jetsTrue = sorted_by_pt(csTrue.inclusive_jets(jetMinPt));

};


void KinematicsJets::CalculateJetKinematics(fastjet::PseudoJet jet){
  // `jet` is already in the head-on frame, since `jetsRec` was filled with head-on frame momenta
  TLorentzVector pjet(jet.px(), jet.py(), jet.pz(), jet.E());
  TVector3 qTjetVect( vecElectron.Px()+pjet.Px(), vecElectron.Py()+pjet.Py(), 0); // (used only in Lorentz invariant calculations)
  qTjet = qTjetVect.Mag();

  zjet = (vecIonBeam.Dot(pjet))/((vecIonBeam).Dot(vecQ));
  pTjet = jet.pt(); // lab frame pT
  mTjet = jet.mt(); // Transverse Mass
  etajet = jet.eta(); // Jet Pseudorapidity
  phijet = jet.phi(); // Jet Phi
  mjet = jet.m(); // Jet Mass
  ejet = jet.e(); // Jet Energy

  jperp.clear();
  zhad_jet.clear();
  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  int constituentPID = 0; // if we only want zh/jperp for pi+, other tracks

  if(constituentPID == 0){
    for(int i = 0; i < constituents.size(); i++){
      TLorentzVector partVec(constituents[i].px(), constituents[i].py(), constituents[i].pz(), constituents[i].E());
      TVector3 jperpVec = Reject(partVec.Vect(),pjet.Vect());
      jperp.push_back(jperpVec.Mag());
      zhad_jet.push_back( (partVec.Vect()).Mag()/((pjet.Vect()).Mag()) );
    }
  }
  else {
    for(int i = 0; i < constituents.size(); i++){
      TLorentzVector partVec(constituents[i].px(), constituents[i].py(), constituents[i].pz(), constituents[i].E());
      std::map<double,int>::iterator it;
      it = jetConstituents.find(partVec.Px());
      if(it != jetConstituents.end()){
        int pidTrack = it->second;
        if(pidTrack == constituentPID){
          TVector3 jperpVec = Reject(partVec.Vect(),pjet.Vect());
          jperp.push_back(jperpVec.Mag());
          zhad_jet.push_back( (partVec.Vect()).Mag()/((pjet.Vect()).Mag()) );
        }
      }
    }
  }
};


void KinematicsJets::CalculateJetResolution(double deltaRCut){
  // Find Matching Truth-level Jet for each Reco Jet
  double minDeltaR = 100000.0;
  double matchIndex = -1;
  for(unsigned int jt=0; jt<jetsTrue.size(); jt++)
    {
      double deltaEta = etajet - jetsTrue[jt].eta();
      double deltaPhi = TVector2::Phi_mpi_pi(phijet - jetsTrue[jt].phi());
      double deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      if(deltaR < minDeltaR)
	{
	  minDeltaR = deltaR;
	  matchIndex = jt;
	}
    }

  // Pass DeltaR to Saved Variables
  if(matchIndex == -1)
    deltaRjet = -1.0;
  else
    deltaRjet = minDeltaR;

  // Calculate Resolution for Matched Jets
  if(minDeltaR < deltaRCut && matchIndex != -1) matchStatusjet = 0; // Matched Jet Found
  if(minDeltaR > deltaRCut && matchIndex != -1) matchStatusjet = 1; // True Jet Found, but Outside Match Criteria
  if(matchIndex == -1) matchStatusjet = 2; // No True Jet Found
  if(minDeltaR < deltaRCut && matchIndex != -1)
    {
      pTmtjet = jetsTrue[matchIndex].pt(); // lab frame matched true jet pT
      mTmtjet = jetsTrue[matchIndex].mt(); // Matched true Jet Transverse Mass
      etamtjet = jetsTrue[matchIndex].eta(); // Matched true Jet Pseudorapidity
      phimtjet = jetsTrue[matchIndex].phi(); // Matched true Jet Phi
      mmtjet = jetsTrue[matchIndex].m(); // Matched True Jet Mass
      emtjet= jetsTrue[matchIndex].e(); // Matched true Jet Energy
    }
};


#ifdef INCLUDE_CENTAURO
void KinematicsJets::GetBreitFrameJets(
  TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
  TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
  )
{
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();
  itParticle.Reset();
  std::vector<fastjet::PseudoJet> particles;
  std::vector<fastjet::PseudoJet> particlesTrue;

  jetConstituents.clear();

  double highPT = -1;
  TLorentzVector eleTrue;
  while(GenParticle* part = (GenParticle*) itParticle()){
    if(part->PID == 11){
      if(part->PT > highPT){
        highPT = part->PT;
        eleTrue = part->P4();
      }
    }
  }
  TLorentzVector vecQTrue = vecEleBeam - eleTrue;
  double Q2true = -1*vecQTrue.M2();
  double xtrue = Q2true / ( 2 * vecQTrue.Dot(vecIonBeam) );

  TLorentzVector breitVecTrue = vecQTrue + 2*xtrue*vecIonBeam;
  TLorentzVector breitVec = vecQ + 2*x*vecIonBeam;
  TVector3 breitBoostTrue = -1*breitVecTrue.BoostVector();
  TVector3 breitBoost = -1*breitVec.BoostVector();
  itParticle.Reset();

  while(Track *eflowTrack = (Track*)itEFlowTrack() ){
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    if(!isnan(eflowTrackp4.E()) && eflowTrackp4 != vecElectron){
      if(std::abs(eflowTrack->Eta) < 4.0 && eflowTrack->PT > 0.2){
        eflowTrackp4.Boost(breitBoost);
        particles.push_back(fastjet::PseudoJet(eflowTrackp4.Px(),eflowTrackp4.Py(),eflowTrackp4.Pz(),eflowTrackp4.E()));

        GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
        TLorentzVector partp4 = trackParticle->P4();
        partp4.Boost(breitBoostTrue);
        particlesTrue.push_back(fastjet::PseudoJet(partp4.Px(),partp4.Py(),partp4.Pz(),partp4.E()));

        jetConstituents.insert({eflowTrackp4.Px(), eflowTrack->PID});

      }
    }
  }
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0 &&
          sqrt(towerPhotonp4.Px()*towerPhotonp4.Px()+towerPhotonp4.Py()*towerPhotonp4.Py()) > 0.2
          )
      {
        towerPhotonp4.Boost(breitBoost);
        particles.push_back(fastjet::PseudoJet(towerPhotonp4.Px(),towerPhotonp4.Py(),towerPhotonp4.Pz(),towerPhotonp4.E()));

        for(int i = 0; i < towerPhoton->Particles.GetEntries(); i++){
          GenParticle *photonPart = (GenParticle*)towerPhoton->Particles.At(i);
          TLorentzVector photonp4 = photonPart->P4();
          photonp4.Boost(breitBoostTrue);
          particlesTrue.push_back(fastjet::PseudoJet(photonp4.Px(),photonp4.Py(),photonp4.Pz(),photonp4.E()));
        }
      }
    }
  }
  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if( !isnan(towerNeutralHadronp4.E()) &&
        sqrt(towerNeutralHadronp4.Px()*towerNeutralHadronp4.Px()+towerNeutralHadronp4.Py()*towerNeutralHadronp4.Py()) > 0.2
        )
    {
      if( std::abs(towerNeutralHadron->Eta) < 4.0 ){
        towerNeutralHadronp4.Boost(breitBoost);
        particles.push_back(
            fastjet::PseudoJet(towerNeutralHadronp4.Px(),towerNeutralHadronp4.Py(),towerNeutralHadronp4.Pz(),towerNeutralHadronp4.E())
          );

        for(int i = 0; i < towerNeutralHadron->Particles.GetEntries(); i++){
          GenParticle *nhadPart = (GenParticle*)towerNeutralHadron->Particles.At(i);
          TLorentzVector nhadp4 = nhadPart->P4();
          nhadp4.Boost(breitBoostTrue);
          particlesTrue.push_back(fastjet::PseudoJet(nhadp4.Px(),nhadp4.Py(),nhadp4.Pz(),nhadp4.E()));
        }
      }
    }
  }

  double R = 0.8;
  contrib::CentauroPlugin centauroPlugin(R);
  fastjet::JetDefinition jet_def(&centauroPlugin);

  csRec = fastjet::ClusterSequence(particles, jet_def);
  csTrue = fastjet::ClusterSequence(particlesTrue, jet_def);
  breitJetsRec = sorted_by_pt(csRec.inclusive_jets());
  breitJetsTrue = sorted_by_pt(csTrue.inclusive_jets());
};


void KinematicsJets::CalculateBreitJetKinematics(fastjet::PseudoJet jet){
  TLorentzVector pjet(jet.px(), jet.py(), jet.pz(), jet.E());
  TLorentzVector pjetLab = pjet;

  TLorentzVector breitVec = vecQ + 2*x*vecIonBeam;
  TVector3 breitBoost = -1*breitVec.BoostVector();

  pjetLab.Boost(-1*breitBoost);
  pTjet = sqrt(pjetLab.Px()*pjetLab.Px() + pjetLab.Py()*pjetLab.Py());

  TLorentzVector vecElectronBreit = vecElectron;
  vecElectronBreit.Boost(breitBoost);
  TVector3 qTjetVect(vecElectronBreit.Px()+pjet.Px(), vecElectronBreit.Py()+pjet.Py(), 0);
  qTjet = qTjetVect.Mag();

  TLorentzVector nbreit(0,0,1/sqrt(Q2),1/sqrt(Q2));
  double zjet = nbreit*pjet;

  jperp.clear();
  zhad_jet.clear();
  std::vector<fastjet::PseudoJet> constituents = jet.constituents();
  int constituentPID = 0; // if we only want zh/jperp for pi+, other tracks

  if(constituentPID == 0){
    for(int i = 0; i < constituents.size(); i++){
      TLorentzVector partVec(constituents[i].px(), constituents[i].py(), constituents[i].pz(), constituents[i].E());
      TVector3 jperpVec = Reject(partVec.Vect(),pjet.Vect());
      jperp.push_back(jperpVec.Mag());
      zhad_jet.push_back( (partVec.Vect()).Mag()/((pjet.Vect()).Mag()) );
    }
  }
  else{
    for(int i = 0; i < constituents.size(); i++){
      TLorentzVector partVec(constituents[i].px(), constituents[i].py(), constituents[i].pz(), constituents[i].E());
      std::map<double,int>::iterator it;
      it = jetConstituents.find(partVec.Px());
      if(it != jetConstituents.end()){
        int pidTrack = it->second;
        if( pidTrack == constituentPID){
          TVector3 jperpVec = Reject(partVec.Vect(),pjet.Vect());
          jperp.push_back(jperpVec.Mag());
          zhad_jet.push_back( (partVec.Vect()).Mag()/((pjet.Vect()).Mag()) );
        }
      }
    }
  }

};
#endif // ifdef INCLUDE_CENTAURO
// end DELPHES-only methods //////////////////////////////////////////////////////


KinematicsJets::~KinematicsJets() {};
#endif // ifndef EXCLUDE_DELPHES
