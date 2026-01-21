// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Connor Pecar, Duane Byer, Sanghwa Park, Brian Page

/* NOTE:
 * if you make changes, MAINTAIN DOCUMENTATION IN ../doc/kinematics.md
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ROOT
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TRandomGen.h>
#include <TClonesArray.h>

// Delphes
#ifndef EXCLUDE_DELPHES
#include <classes/DelphesClasses.h>
#endif
// pybind (for ML models using python packages)
#ifdef SIDIS_MLPRED
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#endif

using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::isnan;

class Kinematics
{
  public:
    Kinematics(Double_t enEleBeam, Double_t enIonBeam, Double_t crossAng);
    ~Kinematics();

    // SIDIS calculators
    Bool_t CalculateDIS(TString recmethod); // return true if succeeded
    void CalculateHadronKinematics();

    // hadronic final state (HFS)
    void AddToHFS(TLorentzVector p4_);
    void AddToHFSTree(TLorentzVector p4, int pid);   	      
    void AddTrackToHFSTree(TLorentzVector p4, int pid);			   
    void SubtractElectronFromHFS();
    void ResetHFS();


    // DELPHES-specific methods //////////////////////////
#ifndef EXCLUDE_DELPHES

    // hadronic final state (HFS)
    void GetHFS(
        TObjArrayIter itTrack,
        TObjArrayIter itEFlowTrack,
        TObjArrayIter itEFlowPhoton,
        TObjArrayIter itEFlowNeutralHadron,
        TObjArrayIter itpfRICHTrack,
        TObjArrayIter itDIRCepidTrack,   TObjArrayIter itDIRChpidTrack,
        TObjArrayIter itBTOFepidTrack,   TObjArrayIter itBTOFhpidTrack,
        TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
        );
    void GetTrueHFS(TObjArrayIter itParticle);

    // PID
    int GetTrackPID(
        Track *track,
        TObjArrayIter itpfRICHTrack,
        TObjArrayIter itDIRCepidTrack, TObjArrayIter itDIRChpidTrack,
        TObjArrayIter itBTOFepidTrack, TObjArrayIter itBTOFhpidTrack,
        TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
        );

#endif // ifndef EXCLUDE_DELPHES
    // end DELPHES-specific methods //////////////////////////


    // kinematics (should be Double_t, if going in SidisTree)
    Double_t W,Q2,Nu,x,y,s; // DIS
    Double_t pLab,pTlab,phiLab,etaLab,z,pT,qT,mX,xF,phiH,phiS; // hadron
    Double_t sigmah, Pxh, Pyh; // hadronic final state variables
    TLorentzVector hadronSumVec;
  
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

    // HFS tree objects
    Int_t nHFS;    
    std::vector<double> hfspx;
    std::vector<double> hfspy;
    std::vector<double> hfspz;
    std::vector<double> hfsE;
    std::vector<int> hfspid;

    // HFS tree select FS tracks and matched true p4
    // (NOT full true HFS)
    std::vector<double> selectedHadPx;
    std::vector<double> selectedHadPy;
    std::vector<double> selectedHadPz;
    std::vector<double> selectedHadE;
    std::vector<int> selectedHadPID;

    // TMVA for ML sidis reconstruction
    std::vector<std::vector<float>> hfsinfo;
    std::vector<float> globalinfo;

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
    void CalculateDISbyML();
    // calculate 4-momenta components of q and W (`vecQ` and `vecW`) as well as
    // derived invariants `W` and `nu`
    void GetQWNu_electronic();
    void GetQWNu_hadronic();
    void GetQWNu_quadratic();
    void GetQWNu_ML();
  private:
    static const Int_t asymInjectN = 2;
    Double_t moduVal[asymInjectN];
    Double_t ampVal[asymInjectN];
    Double_t asymInject;
    std::unique_ptr<TRandom> RNG;
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
#ifdef SIDIS_MLPRED
    py::object keras, tensorflow;
    py::object efnpackage;
    py::function pfnimport;
    py::object model;
    py::object modelload;
    std::string modelname = "pfn_testEpic_000-2_vecQele_nHFS2_500_bs10k_bestValLoss";
#endif  

  ClassDef(Kinematics,1);
};
