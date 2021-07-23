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
  s = (vecEleBeam+vecIonBeam).M2();

  // default transverse spin (needed for phiS calculation)
  tSpin = 1; // +1=spin-up, -1=spin-down

};


// calculates q,W, boost vecs from quadratic formula
void Kinematics::getqWQuadratic(){
  float a = 1 - (vecIonBeam.E()/vecIonBeam.Pz())*(vecIonBeam.E()/vecIonBeam.Pz());
  float b = 2*y*vecIonBeam.E()*(vecIonBeam*vecEleBeam)/(vecIonBeam.Pz()*vecIonBeam.Pz());
  float c = Q2 - Pxh*Pxh - Pyh*Pyh - pow(y*(vecIonBeam*vecEleBeam)/(vecIonBeam.Pz()),2.0);
  
  float qE1, qE2, qE, qz;
  if(b*b>4*a*c && a != 0){
    qE1 = (-1*b+sqrt(b*b-4*a*c))/(2*a);
    qE2 = (-1*b-sqrt(b*b-4*a*c))/(2*a);
    
    if( (qE1) < (qE2) ) qE = qE1;
    else qE = qE2;
    
    qz = (-1*y*(vecEleBeam*vecIonBeam)+vecIonBeam.E()*qE)/(vecIonBeam.Pz());
    
    vecQ.SetPxPyPzE(Pxh, Pyh, qz, qE);
    vecW = vecIonBeam - vecQ;
    W = vecW.M();
    Nu = vecIonBeam.Dot(vecQ)/IonMass;
    this->SetBoostVecs();                                                                                                                                                                                                        
  }   
};


// calculate DIS kinematics using scattered electron
// - needs `vecElectron` set
void Kinematics::CalculateDISbyElectron() {
  vecW = vecEleBeam + vecIonBeam - vecElectron;
  vecQ = vecEleBeam - vecElectron;
  W = vecW.M();
  Q2 = -1 * vecQ.M2();
  Nu = vecIonBeam.Dot(vecQ) / IonMass;
  x = Q2 / ( 2 * vecQ.Dot(vecIonBeam) );
  y = vecIonBeam.Dot(vecQ) / vecIonBeam.Dot(vecEleBeam);
  this->SetBoostVecs();
};

// calculate DIS kinematics using JB method
// sets q, W using quadratic equation
void Kinematics::CalculateDISbyJB(){
  y = sigmah/(2*vecEleBeam.E());
  Q2 = (Pxh*Pxh + Pyh*Pyh)/(1-y);
  x = Q2/(s*y);
  Kinematics::getqWQuadratic();
};

// calculate DIS kinematics using DA method                                                                                                                                                                                                  
// sets q, W using quadratic equation
// requires 'vecElectron' set
void Kinematics::CalculateDISbyDA(){
    float thetah = acos( (Pxh*Pxh+Pyh*Pyh - sigmah*sigmah)/(Pxh*Pxh+Pyh*Pyh+sigmah*sigmah) );
    float thetae = vecElectron.Theta();
    Q2 = 4.0*vecEleBeam.E()*vecEleBeam.E()*sin(thetah)*(1+cos(thetae))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
    y = (sin(thetae)*(1-cos(thetah)))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
    x = Q2/(s*y);
    Kinematics::getqWQuadratic();    
};

// calculate DIS kinematics using mixed method                                                                                                                                                                                    
// requires 'vecElectron' set                                                                                                                                                                                                        

void Kinematics::CalculateDISbyMixed(){
  vecQ = vecEleBeam - vecElectron;
  Q2 = -1*vecQ.M2();
  y = sigmah/(2*vecEleBeam.E());
  x = Q2/(s*y);
  vecW = vecEleBeam + vecIonBeam - vecElectron;
  W = vecW.M();
  Nu = vecIonBeam.Dot(vecQ)/IonMass;
  this->SetBoostVecs();                                                                                                                                                                                                            
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
  // feynman-x
  xF = 2 * CvecHadron.Vect().Dot(CvecQ.Vect()) /
      (W * CvecQ.Vect().Mag());
  // phiH
  phiH = PlaneAngle(
      IvecQ.Vect(), IvecElectron.Vect(),
      IvecQ.Vect(), IvecHadron.Vect()
      );
  // phiS
  vecSpin.SetXYZT(0,tSpin,0,0); // Pauli-Lubanski pseudovector
  //this->BoostToBreitFrame(vecSpin,IvecSpin); // TODO: check if other frames matter
  phiS = PlaneAngle(
      IvecQ.Vect(), IvecElectron.Vect(),
      IvecQ.Vect(), vecSpin.Vect()
      );
  // pT, in perp frame (transverse to q), in ion rest frame
  pT = Reject(
      IvecHadron.Vect(),
      IvecQ.Vect()
      ).Mag();
  // qT
  qT = pT / z;
};

// get PID information from PID systems tracks
int getTrackPID(Track *track, TObjArrayIter itParticle, TObjArrayIter itPIDSystemsTrack){
  GenParticle *trackParticle = (GenParticle*)track->Particle.GetObject();                                                                                                                                            
  GenParticle *detectorParticle;                                                                                                                                                                                                 
  int pidOut = -1;                                                                                                                                                                                                  
  while(Track *detectorTrack = (Track*)itPIDSystemsTrack() ){                                                                                                                                                              
    detectorParticle = (GenParticle*)detectorTrack->Particle.GetObject();                                                                                                                                                             
    if( detectorParticle == trackParticle ) pidOut = detectorTrack->PID;                                                                                                                                               
  }                                                                                                                                                                                                                                         
  return pidOut;
}


// calculates hadronic final state variables from DELPHES tree branches
// expects 'vecElectron' set
void Kinematics::GetHadronicFinalState(TObjArrayIter itTrack, TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton, TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itPIDSystemsTrack, TObjArrayIter itParticle){
  itTrack.Reset();
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();
  itPIDSystemsTrack.Reset();

  itParticle.Reset();
  while(Track *track = (Track*)itTrack() ){  
    TLorentzVector  trackp4 = track->P4();
    if(!isnan(trackp4.E())){
      if( std::abs(track->Eta) >= 4.0  ){ // eta cut?                                                                                                                                                                                    
	int pidTrack = getTrackPID(track, itParticle, itPIDSystemsTrack);
	float trackPt = trackp4.Pt();
	float trackEta = trackp4.Eta();
	float trackPhi = trackp4.Phi();
	float corrmass = correctMass(pidTrack);

	float trackMass = trackp4.M();
	if(corrmass!=0) trackMass = corrmass;
	trackp4.SetPtEtaPhiM(trackPt, trackEta, trackPhi, trackMass);
	
	sigmah += (trackp4.E() - trackp4.Pz());
	Pxh += trackp4.Px();
	Pyh +=trackp4.Py();
      }
    }
  }
  while(Track *eflowTrack = (Track*)itEFlowTrack() ){  
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) < 4.0){
	int pidTrack = getTrackPID(eflowTrack, itParticle, itPIDSystemsTrack);
	float trackPt = eflowTrackp4.Pt();
	float trackEta = eflowTrackp4.Eta();
	float trackPhi = eflowTrackp4.Phi();
	float trackMass = eflowTrackp4.M();
	float corrmass = correctMass(pidTrack);
	if(corrmass!=0) trackMass = corrmass;
	eflowTrackp4.SetPtEtaPhiM(trackPt, trackEta, trackPhi, trackMass);
	
	sigmah += (eflowTrackp4.E() - eflowTrackp4.Pz());
	Pxh += eflowTrackp4.Px();
	Pyh += eflowTrackp4.Py();
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
      }
    }
  }
  if(!isnan(vecElectron.E())){
    sigmah -= (vecElectron.E() - vecElectron.Pz());
    Pxh -= vecElectron.Px();
    Pyh -= vecElectron.Py();
  }           
};


Kinematics::~Kinematics() {
};

