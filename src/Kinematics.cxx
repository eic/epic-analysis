#include "Kinematics.h"

ClassImp(Kinematics)


Kinematics::Kinematics(
    Double_t enEleBeam, /*GeV*/
    Double_t enIonBeam, /*GeV*/
    Double_t crossAng /*mrad*/
    ) {

  // convert crossing angle to rad
  crossAng *= 1e-3;

  // set ion mass
  IonMass = ProtonMass();

  // set beam 4-momenta // TODO: get proper beams from Brian
  Double_t momEleBeam = EMtoP(enEleBeam,ElectronMass());
  Double_t momIonBeam = EMtoP(enIonBeam,IonMass);
  vecEleBeam.SetPxPyPzE(
      0,
      0,
      -momEleBeam,
      enEleBeam
      );
  vecIonBeam.SetPxPyPzE(
      momIonBeam * TMath::Sin(crossAng),
      0,
      momIonBeam * TMath::Cos(crossAng),
      enIonBeam
      );

};


// calculate DIS kinematics using scattered electron
void Kinematics::DISbyElectron() {
  vecW = vecEleBeam + vecIonBeam - vecElectron; // TODO: cross check this with dispin and eiscim
  vecQ = vecEleBeam - vecElectron;
  W = vecW.M();
  Q2 = -1 * vecQ.M2();
  Nu = vecIonBeam.Dot(vecQ) / IonMass;
  x = Q2 / ( 2 * vecQ.Dot(vecIonBeam) );
  y = vecIonBeam.Dot(vecQ) / vecIonBeam.Dot(vecEleBeam);
};


Kinematics::~Kinematics() {
};

