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

    // SIDIS calculators
    Bool_t CalculateDIS(TString recmethod); // return true if succeeded
    void CalculateHadronKinematics();

    // hadronic final state (HFS)
    void GetHFS(
        TObjArrayIter itTrack,
        TObjArrayIter itEFlowTrack,
        TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron,
        TObjArrayIter itpfRICHTrack,
        TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
        );
    void GetTrueHFS(TObjArrayIter itParticle);
    void ResetHFS();
    void SubtractElectronFromHFS();
    void AddToHFS(TLorentzVector p4_);

    // PID
    int getTrackPID(
        Track *track,
	TObjArrayIter itParticle,
        TObjArrayIter itpfRICHTrack,
	TObjArrayIter itbarrelDIRCTrack,
        TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
        );
  
  
    // jet calculators
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
    Double_t sigmah, Pxh, Pyh; // hadronic final state variables
    TLorentzVector hadronSumVec;


    // ADDED BY GREGORY MATOUSEK
    // June 06 2022
    // ---------------------------------------------------------
    Double_t W_e,Q2_e,Nu_e,x_e,y_e,s_e; // Ele DIS
    Double_t W_JB,Q2_JB,Nu_JB,x_JB,y_JB,s_JB; // JB  DIS
    Double_t W_DA,Q2_DA,Nu_DA,x_DA,y_DA,s_DA; // DA  DIS
    
    // Ele DIS params
    Double_t e_Ei,e_Ef,e_th;   
    // DA DIS params
    Double_t thetah,thetae;
   


    // ---------------------------------------------------------


    // nucleon transverse spin; if you set this externally,
    // it must be done before calculating `phiS` (before
    // `CalculateHadronKinematics`)
    Int_t tSpin; // should be +1 or -1
    Int_t lSpin; // should be +1 or -1
    Int_t hadPID;

    // polarization
    Double_t polT;
    Double_t polL;
    Double_t polBeam;

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


    // particle masses
    static Double_t ElectronMass() { return 0.000511; };
    static Double_t ProtonMass()   { return 0.938272; };
    static Double_t KaonMass()     { return 0.493677; };
    static Double_t PionMass()     { return 0.139570; };
    Double_t IonMass;


    // lorentz transformations
    // - boost from Lab frame `Lvec` to photon+ion C.o.m. frame `Cvec`
    void BoostToComFrame(TLorentzVector Lvec, TLorentzVector &Cvec);
    // - boost from Lab frame `Lvec` to Ion rest frame `Ivec`
    void BoostToIonFrame(TLorentzVector Lvec, TLorentzVector &Ivec);
    // - boost from Lab frame `Lvec` to ion+electron Beam c.o.m. frame `Bvec`
    void BoostToBeamComFrame(TLorentzVector Lvec, TLorentzVector &Bvec);
    // - tranform from Lab frame `Lvec` to Head-on frame `Hvec`
    void TransformToHeadOnFrame(TLorentzVector Lvec, TLorentzVector &Hvec);
    // transform from Head-on frame `Hvec` back to Lab frame `Lvec`
    void TransformBackToLabFrame(TLorentzVector Hvec, TLorentzVector &Lvec);


    // misc calculations
    // - convert energy,mass to momentum
    static Double_t EMtoP(Double_t energy, Double_t mass) {
      return TMath::Sqrt( TMath::Power(energy,2) - TMath::Power(mass,2) );
    };
    // - vector projection: returns vA projected onto vB
    static TVector3 Project(TVector3 vA, TVector3 vB) {
      if(fabs(vB.Dot(vB))<0.0001) {
        //cerr << "WARNING: Kinematics::Project to null vector" << endl;
        return TVector3(0,0,0);
      };
      return vB * ( vA.Dot(vB) / ( vB.Dot(vB) ) );
    };
    // - vector rejection: returns vC projected onto plane transverse to vD
    static TVector3 Reject(TVector3 vC, TVector3 vD) {
      if(fabs(vD.Dot(vD))<0.0001) {
        //cerr << "WARNING: Kinematics::Reject to null vector" << endl;
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
        //cerr << "WARNING: Kinematics:PlaneAngle (AxB).D=0" << endl;
        return -10000;
      };
      sgn /= fabs(sgn); // sign of (AxB).D
      Double_t numer = crossAB.Dot(crossCD); // (AxB).(CxD)
      Double_t denom = crossAB.Mag() * crossCD.Mag(); // |AxB|*|CxD|
      if(fabs(denom)<0.00001) {
        //cerr << "WARNING: Kinematics:PlaneAngle denominator=0" << endl;
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

    // tests and validation
    void ValidateHeadOnFrame();

    Long64_t countPIDsmeared,countPIDtrue,countHadrons;


    // settings
    Int_t mainFrame;
    enum mainFrame_enum {fLab, fHeadOn};
    Int_t qComponentsMethod;
    enum qComponentsMethod_enum {qQuadratic, qHadronic, qElectronic};
     
  protected:

    // reconstruction methods
    void CalculateDISbyElectron();
    void CalculateDISbyJB();
    void CalculateDISbyDA();
    void CalculateDISbyMixed();
    void CalculateDISbySigma();
    void CalculateDISbyeSigma();

    // calculate 4-momenta components of q and W (`vecQ` and `vecW`) as well as
    // derived invariants `W` and `nu`
    void GetQWNu_electronic();
    void GetQWNu_hadronic();
    void GetQWNu_quadratic();

  private:
    static const Int_t asymInjectN = 2;
    Double_t moduVal[asymInjectN];
    Double_t ampVal[asymInjectN];
    Double_t asymInject;
    TRandom *RNG;
    Float_t RN;
    Bool_t reconOK;

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
    // - head-on frame
    TLorentzVector HvecEleBeam, HvecIonBeam;
    TLorentzVector HvecElectron, HvecW, HvecQ;
    TLorentzVector HvecHadron;
    // - other intermediate frames (for head-on frame transformation)
    TLorentzVector BvecBoost, OvecBoost;
    TVector3 Bboost, Oboost;
    TLorentzVector BvecEleBeam, BvecIonBeam;
    Double_t rotAboutX, rotAboutY;
    // other
    TLorentzVector vecSpin, IvecSpin;


  ClassDef(Kinematics,1);
};

#endif

