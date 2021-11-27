
#include "Kinematics.h"

ClassImp(Kinematics)

Kinematics::Kinematics(
    Double_t enEleBeam, /*GeV*/
    Double_t enIonBeam, /*GeV*/
    Double_t crossAng /*mrad*/
    )
{

  // set ion mass
  IonMass = ProtonMass();

  // revise crossing angle
  crossAng *= 1e-3; // mrad -> rad
  crossAng = -1*TMath::Abs(crossAng); // take -1*abs(crossAng) to enforce the correct sign

  // set beam 4-momenta
  // - electron beam points toward negative z
  // - ion beam points toward positive z, negative x
  // - crossing angle about y-axis is therefore negative, in right-handed coord. system
  // - north is +x, east is +z
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
  s = (vecEleBeam+vecIonBeam).M2();

  // calculate transformations for head-on frame boost
  // - boost lab frame -> c.o.m. frame of proton and ion Beams
  BvecBoost = vecEleBeam + vecIonBeam;
  Bboost = -1*BvecBoost.BoostVector();
  // - boost c.o.m. frame of beams -> back to a frame with energies (nearly) the Original beam energies
  OvecBoost.SetXYZT( 0.0, 0.0, BvecBoost[2], BvecBoost[3] );
  Oboost = OvecBoost.BoostVector();
  // - boost beams to c.o.m. frame of beams
  this->BoostToBeamComFrame(vecEleBeam,BvecEleBeam);
  this->BoostToBeamComFrame(vecIonBeam,BvecIonBeam);
  // - rotation of beams about y to remove x-components
  rotAboutY = -TMath::ATan2( BvecIonBeam.Px(), BvecIonBeam.Pz() );
  BvecEleBeam.RotateY(rotAboutY);
  BvecIonBeam.RotateY(rotAboutY);
  // - rotation of beams about x to remove y-components
  rotAboutX = TMath::ATan2( BvecIonBeam.Py(), BvecIonBeam.Pz() );

  // default transverse spin (needed for phiS calculation)
  tSpin = 1; // +1=spin-up, -1=spin-down
  lSpin = 1;

  // default proton polarization
  polT = 0.80;
  polL = 0.;
  polBeam = 0.;

  // random number generator (for asymmetry injection
  RNG = new TRandomMixMax(91874); // (TODO: fixed seed?)

  countGood=countBad=0; // DEBUG /////////////////////
};


// ---------------------------------
// calculates 4-momenta components of q and W (`vecQ` and `vecW`) as well as
// derived invariants `W` and `nu`
// - use the quadratic method
// - requires `Q2`, `y`, `Pxh`, `Pyh`
void Kinematics::GetQWNu_quadratic(){

  // head on frame // DEBUG ///////
  // /*
  double f = y*(HvecIonBeam.Dot(HvecEleBeam));
  double pz = HvecIonBeam.Pz();
  double py = HvecIonBeam.Py();
  double px = HvecIonBeam.Px();
  double pE = HvecIonBeam.E();
  double hx = HPxh-px;
  double hy = HPyh-py;
  // */

  // lab frame // DEBUG //////////////////
  /*
  double f = y*(vecIonBeam.Dot(vecEleBeam));
  double pz = vecIonBeam.Pz();
  double py = vecIonBeam.Py();
  double px = vecIonBeam.Px();
  double pE = vecIonBeam.E();
  double hx = Pxh-px;
  double hy = Pyh-py;
  */

  double a = 1.0 - (pE*pE)/(pz*pz);
  double b = (2*pE/(pz*pz))*(px*hx + py*hy + f);
  double c = Q2 - hx*hx - hy*hy - (1/(pz*pz))*pow( (f+px*hx+py*hy) ,2.0);

  double qz1, qz2, qE1, qE2, qE, qz;
  if(b*b>=4*a*c && a != 0){
    qE1 = (-1*b+sqrt(b*b-4*a*c))/(2*a);
    qE2 = (-1*b-sqrt(b*b-4*a*c))/(2*a);
    qz1 = (-1*f + pE*qE1 - px*hx - py*hy)/(pz);
    qz2 = (-1*f + pE*qE2 - px*hx - py*hy)/(pz);

    if(fabs(qE1) < fabs(qE2)){
      qE = qE1;
      qz = qz1;
    }
    else{
      qE = qE2;
      qz = qz2;
    }

    // head on frame // DEBUG //////////////
    // /*
    vecQ.SetPxPyPzE(hx, hy, qz, qE);
    vecW = HvecIonBeam + vecQ;
    W = vecW.M();
    Nu = HvecIonBeam.Dot(vecQ)/IonMass;
    this->TransformBackToLabFrame(vecQ,vecQ); // need to be in lab frame for downstream `CalculateHadronKinematics` (small effect)
    this->TransformBackToLabFrame(vecW,vecW);
    // */

    // lab frame // DEBUG //////////////
    /*
    vecQ.SetPxPyPzE(hx, hy, qz, qE);
    vecW = vecIonBeam + vecQ;
    W = vecW.M();
    Nu = vecIonBeam.Dot(vecQ)/IonMass;
    */

  } else {
    // DEBUG -- this happens a lot more if using lab frame here...
    if(b*b<4*a*c) cerr << "ERROR: negative discriminant in Kinematics::GetQWNu_quadratic; skipping event" << endl;
    else cerr << "ERROR: zero denominator in Kinematics::GetQWNu_quadratic; skipping event" << endl;
    reconOK = false;
  }
};


// calculates 4-momenta components of q and W (`vecQ` and `vecW`) as well as
// derived invariants `W` and `nu`
// - use the hadronic sum 4-momentum `hadronSumVec`
void Kinematics::GetQWNu_hadronic(){
  vecQ = hadronSumVec - vecIonBeam;
  vecW = hadronSumVec;
  W = vecW.M();
  Nu = vecIonBeam.Dot(vecQ)/IonMass;
};


// calculates 4-momenta components of q and W (`vecQ` and `vecW`) as well as
// derived invariants `W` and `nu`
// - use the scattered electron 4-momentum `vecElectron`
void Kinematics::GetQWNu_electronic(){
  vecQ = vecEleBeam - vecElectron;
  vecW = vecEleBeam + vecIonBeam - vecElectron;
  W = vecW.M();
  Nu = vecIonBeam.Dot(vecQ) / IonMass;
};

// ------------------------------------------------------


// function to call different reconstruction methods
Bool_t Kinematics::CalculateDIS(TString recmethod){

  reconOK = true;

  // transform to the head-on frame; not needed by all reconstruction methods,
  // but best to make sure this is done up front
  this->TransformToHeadOnFrame(vecEleBeam,HvecEleBeam);
  this->TransformToHeadOnFrame(vecIonBeam,HvecIonBeam);
  this->TransformToHeadOnFrame(vecElectron,HvecElectron);

  // calculate primary DIS variables, including Q2,x,y,W,nu
  if( recmethod.CompareTo("Ele", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyElectron();
  }
  else if( recmethod.CompareTo("DA", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyDA();
  }
  else if( recmethod.CompareTo("JB", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyJB();
  }
  else if( recmethod.CompareTo("Mixed", TString::kIgnoreCase) == 0){
    this->CalculateDISbyMixed();
  }
  else if( recmethod.CompareTo("Sigma", TString::kIgnoreCase) == 0){
    this->CalculateDISbySigma();
  }
  else if( recmethod.CompareTo("eSigma", TString::kIgnoreCase) == 0){
    this->CalculateDISbyeSigma();
  }
  else {
    cerr << "ERROR: unknown reconstruction method" << endl;
    return false;
  };

  // calculate SIDIS boost vectors
  // - lab frame -> C.o.m. frame of virtual photon and ion
  CvecBoost = vecQ + vecIonBeam;
  Cboost = -1*CvecBoost.BoostVector();
  // - lab frame -> Ion rest frame
  IvecBoost = vecIonBeam;
  Iboost = -1*IvecBoost.BoostVector();

  // calculate depolarization
  // - calculate epsilon, the ratio of longitudinal and transverse photon flux [hep-ph/0611265]
  // - these calculations are Lorentz invariant
  gamma = 2*ProtonMass()*x / TMath::Sqrt(Q2);
  epsilon = ( 1 - y - TMath::Power(gamma*y,2)/4 ) /
    ( 1 - y + y*y/2 + TMath::Power(gamma*y,2)/4 );
  // - factors A,B,C,V,W (see [hep-ph/0611265] using notation from [1408.5721])
  depolA = y*y / (2 - 2*epsilon);
  depolB = depolA * epsilon;
  depolC = depolA * TMath::Sqrt(1-epsilon*epsilon);
  depolV = depolA * TMath::Sqrt(2*epsilon*(1+epsilon));
  depolW = depolA * TMath::Sqrt(2*epsilon*(1-epsilon));
  // - factor ratios (see [1807.10606] eq. 2.3)
  if(depolA==0) depolP1=depolP2=depolP3=depolP4=0;
  else {
    depolP1 = depolB / depolA;
    depolP2 = depolC / depolA;
    depolP3 = depolV / depolA;
    depolP4 = depolW / depolA;
  };

  return reconOK;
};

// calculate DIS kinematics using scattered electron
// - needs `vecElectron` set
void Kinematics::CalculateDISbyElectron() {
  this->GetQWNu_electronic(); // set `vecQ`, `vecW`, `W`, `Nu`
  Q2 = -1 * vecQ.M2();
  x = Q2 / ( 2 * vecQ.Dot(vecIonBeam) );
  y = vecIonBeam.Dot(vecQ) / vecIonBeam.Dot(vecEleBeam);
};

// calculate DIS kinematics using JB method
// sets q, W using quadratic equation
void Kinematics::CalculateDISbyJB(){
  y = Hsigmah/(2*HvecEleBeam.E());
  Q2 = (HPxh*HPxh + HPyh*HPyh)/(1-y);
  x = Q2/(s*y);
  this->GetQWNu_quadratic();
  // this->GetQWNu_hadronic(); // DEBUG -- tends to mess up W ////////////////////////////
  // this->GetQWNu_electronic(); // DEBUG -- does things well ///////////////////////////
};

// calculate DIS kinematics using DA method
// sets q, W using quadratic equation
// requires 'vecElectron' set
void Kinematics::CalculateDISbyDA(){
  // --head on
  float thetah = acos( (HPxh*HPxh+HPyh*HPyh - Hsigmah*Hsigmah)/(HPxh*HPxh+HPyh*HPyh+Hsigmah*Hsigmah) );
  float thetae = HvecElectron.Theta();
  Q2 = 4.0*HvecEleBeam.E()*HvecEleBeam.E()*sin(thetah)*(1+cos(thetae))/(sin(thetah)+sin(thetae)-sin(thetah+thetae)); // DEBUG/////
  y = (sin(thetae)*(1-cos(thetah)))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
  x = Q2/(s*y);
  // --lab // DEBUG ///////////////////////
  // float thetah = acos( (Pxh*Pxh+Pyh*Pyh - sigmah*sigmah)/(Pxh*Pxh+Pyh*Pyh+sigmah*sigmah) );
  // float thetae = vecElectron.Theta();
  // Q2 = 4.0*vecEleBeam.E()*vecEleBeam.E()*sin(thetah)*(1+cos(thetae))/(sin(thetah)+sin(thetae)-sin(thetah+thetae)); // DEBUG/////
  // y = (sin(thetae)*(1-cos(thetah)))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
  // x = Q2/(s*y);

  this->GetQWNu_quadratic();
  // this->GetQWNu_hadronic(); // DEBUG -- tends to mess up W ////////////////////////////
  // this->GetQWNu_electronic(); // DEBUG -- does things well ///////////////////////////
};

// calculate DIS kinematics using mixed method
// requires 'vecElectron' set
void Kinematics::CalculateDISbyMixed(){
  this->GetQWNu_electronic();
  Q2 = -1*vecQ.M2();
  y = Hsigmah/(2*HvecEleBeam.E());
  x = Q2/(s*y);
};

// calculate DIS kinematics using Sigma method
// requires 'vecElectron' set
void Kinematics::CalculateDISbySigma(){
  y = Hsigmah/(Hsigmah + HvecElectron.E()*(1-cos(HvecElectron.Theta())));
  Q2 = (HvecElectron.Px()*HvecElectron.Px() + HvecElectron.Py()*HvecElectron.Py())/(1-y);
  x = Q2/(s*y);
  this->GetQWNu_quadratic();
  // this->GetQWNu_hadronic(); // DEBUG -- tends to mess up W ////////////////////////////
  // this->GetQWNu_electronic(); // DEBUG -- does things well ///////////////////////////
};

// calculate DIS kinematics using eSigma method
// requires 'vecElectron' set
void Kinematics::CalculateDISbyeSigma(){
  this->GetQWNu_electronic();
  Q2 = -1*vecQ.M2();
  double ysigma = Hsigmah/(Hsigmah + HvecElectron.E()*(1-cos(HvecElectron.Theta())));
  double Q2sigma = (HvecElectron.Px()*HvecElectron.Px() + HvecElectron.Py()*HvecElectron.Py())/(1-y);
  double xsigma = Q2sigma/(s*ysigma);    
  y = Q2/(s*xsigma);
  x = xsigma;
  this->GetQWNu_quadratic();
  // this->GetQWNu_hadronic(); // DEBUG -- tends to mess up W ////////////////////////////
  // this->GetQWNu_electronic(); // DEBUG -- does things well ///////////////////////////
};



// calculate hadron kinematics
// - calculate DIS kinematics first, so we have `vecQ`, etc.
// - needs `vecHadron` set
void Kinematics::CalculateHadronKinematics() {
  // hadron momentum
  pLab = vecHadron.P();
  pTlab = vecHadron.Pt();
  phiLab = vecHadron.Phi();
  etaLab = vecHadron.Eta();
  // hadron z
  z = vecIonBeam.Dot(vecHadron) / vecIonBeam.Dot(vecQ);
  // missing mass
  mX = (vecW-vecHadron).M(); // missing mass
  // boosts
  this->BoostToComFrame(vecHadron,CvecHadron);
  this->BoostToComFrame(vecQ,CvecQ);
  this->BoostToIonFrame(vecHadron,IvecHadron);
  this->BoostToIonFrame(vecQ,IvecQ);
  this->BoostToIonFrame(vecElectron,IvecElectron);
  // feynman-x: calculated in photon+ion c.o.m. frame
  xF = 2 * CvecHadron.Vect().Dot(CvecQ.Vect()) /
      (W * CvecQ.Vect().Mag());
  // phiH: calculated in ion rest frame
  phiH = AdjAngle(PlaneAngle(
      IvecQ.Vect(), IvecElectron.Vect(),
      IvecQ.Vect(), IvecHadron.Vect()
      ));
  // phiS: calculated in ion rest frame
  tSpin = RNG->Uniform() < 0.5 ? 1 : -1;
  lSpin = RNG->Uniform() < 0.5 ? 1 : -1;
  vecSpin.SetXYZT(0,1,0,0); // Pauli-Lubanski pseudovector, in lab frame
  this->BoostToIonFrame(vecSpin,IvecSpin); // boost to ion rest frame
  phiS = AdjAngle(PlaneAngle(
      IvecQ.Vect(), IvecElectron.Vect(),
      IvecQ.Vect(), IvecSpin.Vect()
      ));
  // pT, in perp frame (transverse to q): calculated in ion rest frame
  pT = Reject(
      IvecHadron.Vect(),
      IvecQ.Vect()
      ).Mag();
  // qT
  qT = pT / z;
};

// validate transformations to the head-on frame
void Kinematics::ValidateHeadOnFrame() {
  this->BoostToIonFrame(vecEleBeam,IvecEleBeam);
  this->BoostToIonFrame(vecIonBeam,IvecIonBeam);
  this->TransformToHeadOnFrame(vecIonBeam,HvecIonBeam);
  this->TransformToHeadOnFrame(vecEleBeam,HvecEleBeam);
  this->TransformToHeadOnFrame(vecIonBeam,HvecIonBeam);
  this->TransformToHeadOnFrame(vecElectron,HvecElectron);
  this->TransformToHeadOnFrame(vecHadron,HvecHadron);
  printf("\nVALIDATION:\n");
  printf("lab E:     "); vecEleBeam.Print();
  printf("lab I:     "); vecIonBeam.Print();
  printf("ion RF  E: "); IvecEleBeam.Print();
  printf("ion RF  I: "); IvecIonBeam.Print();
  printf("head-on E: "); HvecEleBeam.Print();
  printf("head-on I: "); HvecIonBeam.Print();
  printf("---\n");
  printf("lab electron:     "); vecElectron.Print();
  printf("head-on electron: "); HvecElectron.Print();
  printf("difference:       "); (vecElectron-HvecElectron).Print();
  printf("---\n");
  printf("lab hadron:     "); vecHadron.Print();
  printf("head-on hadron: "); HvecHadron.Print();
  printf("difference:     "); (vecHadron-HvecHadron).Print();
};


// get PID information from PID systems tracks
int Kinematics::getTrackPID(
    Track *track,
    TObjArrayIter itpfRICHTrack,
    TObjArrayIter itDIRCepidTrack, TObjArrayIter itDIRChpidTrack,
    TObjArrayIter itBTOFepidTrack, TObjArrayIter itBTOFhpidTrack,
    TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
    ) {

  itpfRICHTrack.Reset();
  itDIRCepidTrack.Reset();   itDIRChpidTrack.Reset();
  itBTOFepidTrack.Reset();   itBTOFhpidTrack.Reset();
  itdualRICHagTrack.Reset(); itdualRICHcfTrack.Reset();
  GenParticle *trackParticle = (GenParticle*)track->Particle.GetObject();
  GenParticle *detectorParticle;

  // TODO: make this less repetitive:

  while(Track *detectorTrack = (Track*)itpfRICHTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }

  while(Track *detectorTrack = (Track*)itDIRCepidTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }
  while(Track *detectorTrack = (Track*)itDIRChpidTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }

  while(Track *detectorTrack = (Track*)itBTOFepidTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }
  while(Track *detectorTrack = (Track*)itBTOFhpidTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }

  while(Track *detectorTrack = (Track*)itdualRICHagTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }
  while(Track *detectorTrack = (Track*)itdualRICHcfTrack() ){
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();
    if( detectorParticle == trackParticle ) return detectorTrack->PID;
  }

  return -1; // not found
}


// calculates hadronic final state variables from DELPHES tree branches
// expects 'vecElectron' set
// - calculates `sigmah`, `Pxh`, and `Pyh` in the lab and head-on frames
void Kinematics::GetHadronicFinalState(
    TObjArrayIter itTrack,
    TObjArrayIter itEFlowTrack,
    TObjArrayIter itEFlowPhoton,
    TObjArrayIter itEFlowNeutralHadron,
    TObjArrayIter itpfRICHTrack,
    TObjArrayIter itDIRCepidTrack,   TObjArrayIter itDIRChpidTrack,
    TObjArrayIter itBTOFepidTrack,   TObjArrayIter itBTOFhpidTrack,
    TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
    ) {

  itTrack.Reset();
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();

  sigmah = Pxh = Pyh = 0; // lab frame
  Hsigmah = HPxh = HPyh = 0; // head-on frame

  hadronSumVec.SetPxPyPzE(0,0,0,0); // lab-frame

  while(Track *track = (Track*)itTrack() ){
    TLorentzVector  trackp4 = track->P4();
    if(!isnan(trackp4.E())){
      if( std::abs(track->Eta) < 4.0  ){

        int pid = getTrackPID( // get smeared PID
            track,
            itpfRICHTrack,
            itDIRCepidTrack, itDIRChpidTrack,
            itBTOFepidTrack, itBTOFhpidTrack,
            itdualRICHagTrack, itdualRICHcfTrack
            );
        if(pid != -1){ // if PID determined, set mass of `trackp4` accordingly
          trackp4.SetPtEtaPhiM(trackp4.Pt(),trackp4.Eta(),trackp4.Phi(),correctMass(pid));
          countGood++; // DEBUG /////////////////////
        }
        else { // otherwise, true PID is assumed
          countBad++; // DEBUG /////////////////////
          //continue; // DEBUG /////////////////////
        }

        sigmah += (trackp4.E() - trackp4.Pz());
        Pxh += trackp4.Px();
        Pyh += trackp4.Py();
        hadronSumVec += trackp4;
        this->TransformToHeadOnFrame(trackp4,trackp4);
        Hsigmah += (trackp4.E() - trackp4.Pz());
        HPxh += trackp4.Px();
        HPyh += trackp4.Py();
      }
    }    
  }

  while(Track *eflowTrack = (Track*)itEFlowTrack() ){
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) >= 4.0){
        sigmah += (eflowTrackp4.E() - eflowTrackp4.Pz());
        Pxh += eflowTrackp4.Px();
        Pyh += eflowTrackp4.Py();
        hadronSumVec += eflowTrackp4;
        this->TransformToHeadOnFrame(eflowTrackp4,eflowTrackp4);
        Hsigmah += (eflowTrackp4.E() - eflowTrackp4.Pz());
        HPxh += eflowTrackp4.Px();
        HPyh += eflowTrackp4.Py();
      }
    }
  }
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0  ){
        sigmah += (towerPhotonp4.E() - towerPhotonp4.Pz());
        Pxh += towerPhotonp4.Px();
        Pyh += towerPhotonp4.Py();
        hadronSumVec += towerPhotonp4;
        this->TransformToHeadOnFrame(towerPhotonp4,towerPhotonp4);
        Hsigmah += (towerPhotonp4.E() - towerPhotonp4.Pz());
        HPxh += towerPhotonp4.Px();
        HPyh += towerPhotonp4.Py();
      }
    }
  }
  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if(!isnan(towerNeutralHadronp4.E())){
      if( std::abs(towerNeutralHadron->Eta) < 4.0 ){
        sigmah += (towerNeutralHadronp4.E() - towerNeutralHadronp4.Pz());
        Pxh += towerNeutralHadronp4.Px();
        Pyh += towerNeutralHadronp4.Py();
        hadronSumVec += towerNeutralHadronp4;
        this->TransformToHeadOnFrame(towerNeutralHadronp4,towerNeutralHadronp4);
        Hsigmah += (towerNeutralHadronp4.E() - towerNeutralHadronp4.Pz());
        HPxh += towerNeutralHadronp4.Px();
        HPyh += towerNeutralHadronp4.Py();
      }
    }
  }
  
  if(!isnan(vecElectron.E())){
    sigmah -= (vecElectron.E() - vecElectron.Pz());
    Pxh -= vecElectron.Px();
    Pyh -= vecElectron.Py();
    hadronSumVec -= vecElectron;
    this->TransformToHeadOnFrame(vecElectron,HvecElectron);
    Hsigmah -= (HvecElectron.E() - HvecElectron.Pz());
    HPxh -= HvecElectron.Px();
    HPyh -= HvecElectron.Py();
  }
};

void Kinematics::GetHadronicFinalStateTrue(TObjArrayIter itParticle){
  itParticle.Reset();

  Hsigmah = 0;
  HPxh = 0;
  HPyh = 0;
  sigmah = 0;
  Pxh = 0;
  Pyh = 0;
  hadronSumVec.SetPxPyPzE(0,0,0,0);

  while(GenParticle *partTrue = (GenParticle*)itParticle() ){
    if(partTrue->Status == 1){
      TLorentzVector partp4 = partTrue->P4();
      sigmah += (partp4.E() - partp4.Pz());
      Pxh += partp4.Px();
      Pyh += partp4.Py();
      hadronSumVec += partp4;
      this->TransformToHeadOnFrame(partp4,partp4);
      Hsigmah += (partp4.E() - partp4.Pz());
      HPxh += partp4.Px();
      HPyh += partp4.Py();
    }    
  }

  sigmah -= (vecElectron.E()-vecElectron.Pz());
  Pxh -= vecElectron.Px();
  Pyh -= vecElectron.Py();
  hadronSumVec -= vecElectron;
  this->TransformToHeadOnFrame(vecElectron,HvecElectron);
  Hsigmah -= (HvecElectron.E()-HvecElectron.Pz());
  HPxh -= HvecElectron.Px();
  HPyh -= HvecElectron.Py();
};


void Kinematics::GetJets(
    TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton,
    TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle
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
    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) < 4.0 && eflowTrack->PT > 0.1){
        this->TransformToHeadOnFrame(eflowTrackp4,eflowTrackp4);
        particles.push_back(fastjet::PseudoJet(eflowTrackp4.Px(),eflowTrackp4.Py(),eflowTrackp4.Pz(),eflowTrackp4.E()));

        GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
        TLorentzVector partp4 = trackParticle->P4();
        this->TransformToHeadOnFrame(partp4,partp4);
        particlesTrue.push_back(fastjet::PseudoJet(partp4.Px(),partp4.Py(),partp4.Pz(),partp4.E()));

        jetConstituents.insert(std::pair<double,int>(eflowTrackp4.Px(), eflowTrack->PID) );
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

  //double R = 0.8*(M_PI/2.0);
  double R = 0.8;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  csRec = fastjet::ClusterSequence(particles, jet_def);
  csTrue = fastjet::ClusterSequence(particlesTrue, jet_def);
  jetsRec = sorted_by_pt(csRec.inclusive_jets());
  jetsTrue = sorted_by_pt(csTrue.inclusive_jets());

};


#if INCCENTAURO == 1
void Kinematics::GetBreitFrameJets(
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

        jetConstituents.insert(std::pair<double,int>(eflowTrackp4.Px(), eflowTrack->PID) );

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


void Kinematics::CalculateBreitJetKinematics(fastjet::PseudoJet jet){
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
#endif


void Kinematics::CalculateJetKinematics(fastjet::PseudoJet jet){
  // `jet` is already in the head-on frame, since `jetsRec` was filled with head-on frame momenta
  TLorentzVector pjet(jet.px(), jet.py(), jet.pz(), jet.E());
  TVector3 qTjetVect( vecElectron.Px()+pjet.Px(), vecElectron.Py()+pjet.Py(), 0); // (used only in Lorentz invariant calculations)
  qTjet = qTjetVect.Mag();

  zjet = (vecIonBeam.Dot(pjet))/((vecIonBeam).Dot(vecQ));
  pTjet = jet.pt(); // lab frame pT

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


// BOOSTS
/////////////////

// boost from Lab frame `Lvec` to photon+ion C.o.m. frame `Cvec`
void Kinematics::BoostToComFrame(TLorentzVector Lvec, TLorentzVector &Cvec) {
  Cvec=Lvec;
  Cvec.Boost(Cboost);
};

// boost from Lab frame `Lvec` to Ion rest frame `Ivec`
void Kinematics::BoostToIonFrame(TLorentzVector Lvec, TLorentzVector &Ivec) {
  Ivec=Lvec;
  Ivec.Boost(Iboost);
};

// boost from Lab frame `Lvec` to ion+electron Beam c.o.m. frame `Bvec`
void Kinematics::BoostToBeamComFrame(TLorentzVector Lvec, TLorentzVector &Bvec) {
  Bvec=Lvec;
  Bvec.Boost(Bboost);
};

// transform from Lab frame `Lvec` to Head-on frame `Hvec`
void Kinematics::TransformToHeadOnFrame(TLorentzVector Lvec, TLorentzVector &Hvec) {
  this->BoostToBeamComFrame(Lvec,Hvec); // boost to c.o.m. frame of beams
  Hvec.RotateY(rotAboutY); // remove x-component of beams
  Hvec.RotateX(rotAboutX); // remove y-component of beams
  Hvec.Boost(Oboost); // return to frame where beam energies are (nearly) the original
};

// transform from Head-on frame `Hvec` back to Lab frame `Lvec`
void Kinematics::TransformBackToLabFrame(TLorentzVector Hvec, TLorentzVector &Lvec) {
  Lvec=Hvec;
  Lvec.Boost(-1*Oboost); // revert boosts and rotations
  Lvec.RotateX(-rotAboutX);
  Lvec.RotateY(-rotAboutY);
  Lvec.Boost(-1*Bboost); // boost from c.o.m. frame of beams back to lab frame
};


// test a fake asymmetry, for fit code validation
// - assigns `tSpin` based on desired fake asymmetry
void Kinematics::InjectFakeAsymmetry() {
  // modulations, including depolarization factors [1807.10606 eq. 2.2-3]
  moduVal[0] = TMath::Sin(phiH-phiS); // sivers
  moduVal[1] = TMath::Sin(phiH+phiS) * depolP1; // transversity*collins
  // fake amplitudes
  ampVal[0] = 0.1;
  ampVal[1] = 0.1;
  // fake dependence on SIDIS kinematics (linear in x)
  asymInject = 0;
  asymInject +=  ampVal[0]/0.2 * x * moduVal[0];
  asymInject += -ampVal[1]/0.2 * x * moduVal[1];
  //asymInject = ampVal[0]*moduVal[0] + ampVal[1]*moduVal[1]; // constant
  // apply polarization
  asymInject *= polT;
  // generate random number in [0,1]
  RN = RNG->Uniform();
  tSpin = (RN<0.5*(1+asymInject)) ? 1 : -1;
};


Kinematics::~Kinematics() {
};

