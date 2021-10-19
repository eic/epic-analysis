#ifndef Kinematics_
#define Kinematics_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRandomGen.h"

// Delphes
#include "classes/DelphesClasses.h"

// Fastjet
#include "fastjet/ClusterSequence.hh"
#if INCCENTAURO == 1
#include "fastjet/plugins/Centauro/Centauro.hh"
#endif
//using namespace fastjet;

using std::map;
using std::cout;
using std::cerr;
using std::endl;

class Kinematics : public TObject
{
  public:
    Kinematics(Double_t enEleBeam, Double_t enIonBeam, Double_t crossAng);
    ~Kinematics();

    // calculators
    void CalculateDIS(TString recmethod);
    void CalculateDISbyElectron();
    void CalculateDISbyJB();
    void CalculateDISbyDA();
    void CalculateDISbyMixed();
    void CalculateDISbySigma();
    void CalculateDISbyeSigma();
    void getqWQuadratic();
    void CalculateHadronKinematics();
    void GetHadronicFinalState(
        TObjArrayIter itTrack, TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
        );
    void GetJets(
        TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
        );
    void CalculateJetKinematics(fastjet::PseudoJet jet);

    #if INCCENTAURO == 1
    void GetBreitFrameJets(
        TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
        );
    void CalculateBreitJetKinematics(fastjet::PseudoJet jet);
    #endif

    // kinematics (should be Double_t, if going in SimpleTree)
    Double_t W,Q2,Nu,x,y,s; // DIS
    Double_t pLab,pTlab,phiLab,etaLab,z,pT,qT,mX,xF,phiH,phiS; // hadron
    Double_t sigmah, Pxh, Pyh; // hadronic final state

    // nucleon transverse spin; if you set this externally,
    // it must be done before calculating `phiS` (before
    // `CalculateHadronKinematics`)
    Int_t tSpin; // should be +1 or -1

    // polarization
    Double_t pol;

    // depolarization
    Double_t gamma,epsilon;
    // - factors A,B,C,V,W from [hep-ph/0611265] using notation from [1408.5721]
    Double_t depolA, depolB, depolC, depolV, depolW;
    // - ratios of factors, following notation of [1807.10606] eq. 2.3 (cf. eqs. 2.2a,b)
    Double_t depolP1; // for A_UT*sin(phiH+phiS) (collins), A_UT*sin(3phiH-phiS) (pretzelosity)
    Double_t depolP2; // for A_LL*const
    Double_t depolP3; // for twist-3 A_UT
    Double_t depolP4; // for A_LL*cos(phiH)

    // 4-vectors
    // - lab frame
    TLorentzVector vecEleBeam, vecIonBeam;
    TLorentzVector vecElectron, vecW, vecQ;
    TLorentzVector vecHadron;
    // jets
    std::vector<fastjet::PseudoJet> jetsRec, jetsTrue;
    std::vector<fastjet::PseudoJet> breitJetsRec, breitJetsTrue;
    std::map<double, int> jetConstituents;

    fastjet::ClusterSequence csRec;
    fastjet::ClusterSequence csTrue;

    Double_t zjet, pTjet, qTjet;
    std::vector<double> jperp;
    std::vector<double> zhad_jet;
    // struck quark information
    Double_t quarkpT;

    // - c.o.m. frame of virtual photon and ion
    TLorentzVector CvecBoost;
    TVector3 Cboost;
    TLorentzVector CvecEleBeam, CvecIonBeam;
    TLorentzVector CvecElectron, CvecW, CvecQ;
    TLorentzVector CvecHadron;
    // - ion rest frame
    TLorentzVector IvecBoost;
    TVector3 Iboost;
    TLorentzVector IvecEleBeam, IvecIonBeam;
    TLorentzVector IvecElectron, IvecW, IvecQ;
    TLorentzVector IvecHadron;
    // other
    TLorentzVector vecSpin, IvecSpin;


    // particle masses
    static Double_t ElectronMass() { return 0.000511; };
    static Double_t ProtonMass()   { return 0.938272; };
    static Double_t KaonMass()     { return 0.493677; };
    static Double_t PionMass()     { return 0.139570; };
    Double_t IonMass;


    // boost calculations
    // - boost from Lab frame to photon+ion C.o.m. frame
    void BoostToComFrame(TLorentzVector Lvec, TLorentzVector &Cvec) {
      Cvec=Lvec; Cvec.Boost(Cboost); };
    // - boost from Lab frame to Ion rest frame
    void BoostToIonFrame(TLorentzVector Lvec, TLorentzVector &Ivec) {
      Ivec=Lvec; Ivec.Boost(Iboost); };


    // misc calculations
    // - convert energy,mass to momentum
    static Double_t EMtoP(Double_t energy, Double_t mass) {
      return TMath::Sqrt( TMath::Power(energy,2) - TMath::Power(mass,2) );
    };
    // - vector projection: returns vA projected onto vB
    static TVector3 Project(TVector3 vA, TVector3 vB) {
      if(fabs(vB.Dot(vB))<0.0001) {
        cerr << "WARNING: Kinematics::Project to null vector" << endl;
        return TVector3(0,0,0);
      };
      return vB * ( vA.Dot(vB) / ( vB.Dot(vB) ) );
    };
    // - vector rejection: returns vC projected onto plane transverse to vD
    static TVector3 Reject(TVector3 vC, TVector3 vD) {
      if(fabs(vD.Dot(vD))<0.0001) {
        cerr << "WARNING: Kinematics::Reject to null vector" << endl;
        return TVector3(0,0,0);
      };
      return vC - Project(vC,vD);
    };
    // - calculate angle between two planes, spanned by vectors
    static Double_t PlaneAngle(TVector3 vA, TVector3 vB, TVector3 vC, TVector3 vD) {
      TVector3 crossAB = vA.Cross(vB); // AxB
      TVector3 crossCD = vC.Cross(vD); // CxD
      Double_t sgn = crossAB.Dot(vD); // (AxB).D
      if(fabs(sgn)<0.00001) {
        cerr << "WARNING: Kinematics:PlaneAngle (AxB).D=0" << endl;
        return -10000;
      };
      sgn /= fabs(sgn); // sign of (AxB).D
      Double_t numer = crossAB.Dot(crossCD); // (AxB).(CxD)
      Double_t denom = crossAB.Mag() * crossCD.Mag(); // |AxB|*|CxD|
      if(fabs(denom)<0.00001) {
        cerr << "WARNING: Kinematics:PlaneAngle denominator=0" << endl;
        return -10000;
      };
      return sgn * TMath::ACos(numer/denom);
    };
    // - shift angle to the range [-PI,+PI]
    static Double_t AdjAngle(Double_t ang) {
      while(ang>TMath::Pi()) ang-=2*TMath::Pi();
      while(ang<-TMath::Pi()) ang+=2*TMath::Pi();
      return ang;
    };

    // misc. functions for hadronic final state
    float correctMass(int pid){
      float massOut = 0;
      switch(std::abs(pid)){
        case 11:
          massOut = ElectronMass();
          break;
        case 2212:
          massOut = ProtonMass();
          break;
        case 321:
          massOut = KaonMass();
          break;
        case 211:
          massOut = PionMass();
          break;
        default:
          cerr << "ERROR: unknown pid in Kinematics::correctMass" << endl;
          massOut = -1;
      };
      return massOut;
    };

    // asymmetry injection
    void InjectFakeAsymmetry(); // test your own asymmetry, for fit code validation


  private:
    static const Int_t asymInjectN = 2;
    Double_t moduVal[asymInjectN];
    Double_t ampVal[asymInjectN];
    Double_t asymInject;
    TRandom *RNG;
    Float_t RN;

  ClassDef(Kinematics,1);
};

#endif

