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

#include "classes/DelphesClasses.h"


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
    void CalculateDISbyElectron();
    void CalculateDISbyJB();
    void CalculateDISbyDA();
    void CalculateDISbyMixed();
  
    void getqWQuadratic();
    void CalculateHadronKinematics();
    void GetHadronicFinalState(TObjArrayIter itTrack, TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton, TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itPIDSystemsTrack, TObjArrayIter itParticle);
  
    // kinematics
    Double_t W,Q2,Nu,x,y,s; // DIS
    Double_t z,pT,qT,mX,xF,phiH,phiS; // hadron
    Double_t sigmah, Pxh, Pyh; // hadronic final state

    // nucleon transverse spin; if you set this externally,
    // it must be done before calculating `phiS` (before
    // `CalculateHadronKinematics`)
    Int_t tSpin; 

    // 4-vectors
    // - lab frame
    TLorentzVector vecEleBeam, vecIonBeam;
    TLorentzVector vecElectron, vecW, vecQ;
    TLorentzVector vecHadron;
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
    static Double_t KaonMass()   { return 0.493677; };
    static Double_t PionMass()   { return 0.139570; };
    Double_t IonMass;


    // boost calculations
    // - calculate boost vectors
    void SetBoostVecs() {
      // c.o.m. frame of virtual photon and ion
      CvecBoost = vecQ + vecIonBeam;
      Cboost = -1*CvecBoost.BoostVector();
      // ion rest frame
      IvecBoost = vecIonBeam;
      Iboost = -1*IvecBoost.BoostVector();
    };
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
  
    // misc. functions for hadronic final state
    float correctMass(int pid){
      float massOut = 0;                                                                                                                                                                                                             
      switch(std::abs(pid)){
      case 11:
	massOut = ElectronMass();
      case 2212:
	massOut = ProtonMass();
      case 321:
	massOut = KaonMass();
      case 211:
	massOut = PionMass();
      }
      return massOut;
    }


    // CUTS =====================================================
    const Double_t xMinGlobal = 0.05; // minimum x for "large-x"
    Bool_t CutDIS() {
      return x>xMinGlobal /* large x region */
          && W>3.0 /* inelastic region */
          && y>0.00 && y<0.95 /* ymin cut applied elsewhere */
          ;
    };
    Bool_t CutHadron() {
      return z>0.2 && z<0.9
          && vecHadron.Pt()>0.1 /* tracking limit on pT_lab */
          && xF>0 /* bias toward current fragmentation */
          ;
    };
    Bool_t CutFull() {
      return this->CutDIS()
          && this->CutHadron()
          ;
    };
    // ==========================================================


  private:

  ClassDef(Kinematics,1);
};

#endif

