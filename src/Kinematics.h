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
    void DISbyElectron();

    // DIS kinematics
    Double_t W,Q2,Nu,x,y;

    // 4-vectors
    TLorentzVector vecEleBeam, vecIonBeam;
    TLorentzVector vecElectron, vecW, vecQ;

    // misc calculations
    static Float_t EMtoP(Float_t energy, Float_t mass) {
      return TMath::Sqrt( TMath::Power(energy,2) - TMath::Power(mass,2) );
    };

    // particle masses
    static Double_t ElectronMass() { return 0.000511; };
    static Double_t ProtonMass()   { return 0.938272; };
    Double_t IonMass;

  private:

  ClassDef(Kinematics,1);
};

#endif
