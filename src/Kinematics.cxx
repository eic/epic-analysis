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
  double f = y*(vecIonBeam.Dot(vecEleBeam));
  double hx = Pxh;
  double hy = Pyh;
  double pz = vecIonBeam.Pz();
  double py = vecIonBeam.Py();
  double px = vecIonBeam.Px();
  double pE = vecIonBeam.E();
  
  double a = 1.0 - (pE*pE)/(pz*pz);
  double b = (2*pE/(pz*pz))*(px*hx + py*hy + f);
  double c = Q2 - hx*hx - hy*hy - (1/(pz*pz))*pow( (f+px*hx+py*hy) ,2.0);
  
  double qz1, qz2, qE1, qE2, qE, qz;
  if(b*b>4*a*c && a != 0){
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
    
    vecQ.SetPxPyPzE(Pxh, Pyh, qz, qE);
    vecW = vecIonBeam - vecQ;
    W = vecW.M();
    Nu = vecIonBeam.Dot(vecQ)/IonMass;
    this->SetBoostVecs();                                                                                                                                                                                                        
  }   
};

// function to call different reconstruction methods
void Kinematics::CalculateDIS(TString recmethod){
  if( recmethod.CompareTo("Ele", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyElectron();
  }
  if( recmethod.CompareTo("DA", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyDA();
  }
  if( recmethod.CompareTo("JB", TString::kIgnoreCase) == 0 ){
    this->CalculateDISbyJB();
  }
  if( recmethod.CompareTo("Mixed", TString::kIgnoreCase) == 0){
    this->CalculateDISbyMixed();
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
  this->getqWQuadratic();
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
    this->getqWQuadratic();    
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
      if( std::abs(track->Eta) >= 4.0  ){ 
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

void Kinematics::GetJets(TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton, TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle){
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();
  itParticle.Reset();
  
  while(GenParticle *partTrue = (GenParticle*)itParticle() ){
    if( (partTrue->PID == 1 || partTrue->PID == 2) && (partTrue->Status == 23) ){
      // Status: 23->outgoing, but there's also 63->outgoing beam remnant. Which do we want?
      // from pythia 8 documentation
      quarkpT = partTrue->PT;
    }    
  }
  
  std::vector<PseudoJet> particles;
  std::vector<PseudoJet> particlesTrue;
  jetConstituents.clear();
  // looping over final state particles, adding to particles vector 
  while(Track *eflowTrack = (Track*)itEFlowTrack() ){
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) < 4.0 && eflowTrack->PT > 0.1){
	particles.push_back(PseudoJet(eflowTrackp4.Px(),eflowTrackp4.Py(),eflowTrackp4.Pz(),eflowTrackp4.E()));

	GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
	TLorentzVector partp4 = trackParticle->P4();	
	particlesTrue.push_back(PseudoJet(partp4.Px(),partp4.Py(),partp4.Pz(),partp4.E()));

	jetConstituents.insert(std::pair<double,int>(eflowTrackp4.Px(), eflowTrack->PID) );
      }
    }
  }
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0){
	particles.push_back(PseudoJet(towerPhotonp4.Px(),towerPhotonp4.Py(),towerPhotonp4.Pz(),towerPhotonp4.E()));
	
	for(int i = 0; i < towerPhoton->Particles.GetEntries(); i++){
	  GenParticle *photonPart = (GenParticle*)towerPhoton->Particles.At(i);
	  TLorentzVector photonp4 = photonPart->P4();
	  particlesTrue.push_back(PseudoJet(photonp4.Px(),photonp4.Py(),photonp4.Pz(),photonp4.E()));	  
	}
      }
    }
  }
  
  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if(!isnan(towerNeutralHadronp4.E())){
      if( std::abs(towerNeutralHadron->Eta) < 4.0){
	particles.push_back(PseudoJet(towerNeutralHadronp4.Px(),towerNeutralHadronp4.Py(),towerNeutralHadronp4.Pz(),towerNeutralHadronp4.E()));

        for(int i = 0; i < towerNeutralHadron->Particles.GetEntries(); i++){
          GenParticle *nhadPart = (GenParticle*)towerNeutralHadron->Particles.At(i);
          TLorentzVector nhadp4 = nhadPart->P4();
          particlesTrue.push_back(PseudoJet(nhadp4.Px(),nhadp4.Py(),nhadp4.Pz(),nhadp4.E()));
	}	
      }
    }
  }
  
  //double R = 0.8*(M_PI/2.0);
  double R = 0.8;
  JetDefinition jet_def(antikt_algorithm, R);

  csRec = new ClusterSequence(particles, jet_def);
  csTrue = new ClusterSequence(particlesTrue, jet_def);
  jetsRec = sorted_by_pt(csRec->inclusive_jets());
  jetsTrue = sorted_by_pt(csTrue->inclusive_jets());  

};

#ifdef __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
void Kinematics::GetBreitFrameJets(TObjArrayIter itEFlowTrack, TObjArrayIter itEFlowPhoton, TObjArrayIter itEFlowNeutralHadron, TObjArrayIter itParticle){
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();
  itParticle.Reset();
  std::vector<PseudoJet> particles;
  std::vector<PseudoJet> particlesTrue;

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
        particles.push_back(PseudoJet(eflowTrackp4.Px(),eflowTrackp4.Py(),eflowTrackp4.Pz(),eflowTrackp4.E()));

        GenParticle *trackParticle = (GenParticle*)eflowTrack->Particle.GetObject();
        TLorentzVector partp4 = trackParticle->P4();
	partp4.Boost(breitBoostTrue);
        particlesTrue.push_back(PseudoJet(partp4.Px(),partp4.Py(),partp4.Pz(),partp4.E()));

	jetConstituents.insert(std::pair<double,int>(eflowTrackp4.Px(), eflowTrack->PID) );

      }
    }
  }
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0 && sqrt(towerPhotonp4.Px()*towerPhotonp4.Px()+towerPhotonp4.Py()*towerPhotonp4.Py()) > 0.2){
	towerPhotonp4.Boost(breitBoost);
	particles.push_back(PseudoJet(towerPhotonp4.Px(),towerPhotonp4.Py(),towerPhotonp4.Pz(),towerPhotonp4.E()));
	
        for(int i = 0; i < towerPhoton->Particles.GetEntries(); i++){
          GenParticle *photonPart = (GenParticle*)towerPhoton->Particles.At(i);
          TLorentzVector photonp4 = photonPart->P4();
	  photonp4.Boost(breitBoostTrue);
          particlesTrue.push_back(PseudoJet(photonp4.Px(),photonp4.Py(),photonp4.Pz(),photonp4.E()));
        }
      }
    }
  }
  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if(!isnan(towerNeutralHadronp4.E()) && sqrt(towerNeutralHadronp4.Px()*towerNeutralHadronp4.Px()+towerNeutralHadronp4.Py()*towerNeutralHadronp4.Py()) > 0.2){
      if( std::abs(towerNeutralHadron->Eta) < 4.0 ){
	towerNeutralHadronp4.Boost(breitBoost);
        particles.push_back(PseudoJet(towerNeutralHadronp4.Px(),towerNeutralHadronp4.Py(),towerNeutralHadronp4.Pz(),towerNeutralHadronp4.E()));

        for(int i = 0; i < towerNeutralHadron->Particles.GetEntries(); i++){
          GenParticle *nhadPart = (GenParticle*)towerNeutralHadron->Particles.At(i);
          TLorentzVector nhadp4 = nhadPart->P4();
	  nhadp4.Boost(breitBoostTrue);
          particlesTrue.push_back(PseudoJet(nhadp4.Px(),nhadp4.Py(),nhadp4.Pz(),nhadp4.E()));
        }
      }
    }
  }

  double R = 0.8;
  contrib::CentauroPlugin centauroPlugin(R);
  JetDefinition jet_def(&centauroPlugin);
  
  csRec = new ClusterSequence(particles, jet_def);
  csTrue = new ClusterSequence(particlesTrue, jet_def);
  breitJetsRec = sorted_by_pt(csRec->inclusive_jets());
  breitJetsTrue = sorted_by_pt(csTrue->inclusive_jets());
  
};

void Kinematics::CalculateBreitJetKinematics(PseudoJet jet){
  TLorentzVector pjet(jet.px(), jet.py(), jet.pz(), jet.E());
  TLorentzVector pjetLab = pjet;

  TLorentzVector breitVec = vecQ + 2*x*vecIonBeam;
  TVector3 breitBoost = -1*breitVec.BoostVector();

  pjetLab.Boost(-1*breitBoost);
  pTjet = sqrt(pjetLab.Px()*pjetLab.Px() + pjetLab.Py()*pjetLab.Py());

  TLorentzVector vecElectronBreit = vecElectron;
  vecElectronBreit.Boost(breitBoost);
  TVector3 qT(vecElectronBreit.Px()+pjet.Px(), vecElectronBreit.Py()+pjet.Py(), 0);
  qTjet = qT.Mag();

  TLorentzVector nbreit(0,0,1/sqrt(Q2),1/sqrt(Q2));
  double zjet = nbreit*pjet;

  jperp.clear();
  zhad_jet.clear();
  std::vector<PseudoJet> constituents = jet.constituents();
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

void Kinematics::CalculateJetKinematics(PseudoJet jet){
  TLorentzVector pjet(jet.px(), jet.py(), jet.pz(), jet.E());
  TVector3 qT( vecElectron.Px()+pjet.Px(), vecElectron.Py()+pjet.Py(), 0);
  qTjet = qT.Mag();

  zjet = (vecIonBeam.Dot(pjet))/((vecIonBeam).Dot(vecQ));
  pTjet = jet.pt(); // lab frame pT

  jperp.clear();
  zhad_jet.clear();
  std::vector<PseudoJet> constituents = jet.constituents();
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
	if(pidTrack == constituentPID){
	  TVector3 jperpVec = Reject(partVec.Vect(),pjet.Vect());
	  jperp.push_back(jperpVec.Mag());
	  zhad_jet.push_back( (partVec.Vect()).Mag()/((pjet.Vect()).Mag()) );
	}
      }
    }
  }
  
};

Kinematics::~Kinematics() {
};

