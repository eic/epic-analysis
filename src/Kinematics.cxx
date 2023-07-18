// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Connor Pecar, Duane Byer, Sanghwa Park, Brian Page

/* NOTE:
 * if you make changes, MAINTAIN DOCUMENTATION IN ../doc/kinematics.md
 */

#include "Kinematics.h"
ClassImp(Kinematics)

Kinematics::Kinematics(
    Double_t enEleBeam, /*GeV*/
    Double_t enIonBeam, /*GeV*/
    Double_t crossAng /*mrad*/    
    )
{
  srand(time(NULL));
  // importing from local python script for ML predictions
  // requires tensorflow, energyflow packages installed
#ifdef SIDIS_MLPRED
  efnpackage = py::module_::import("testEFlowimport");
  pfnimport = efnpackage.attr("eflowPredict");
#endif
  // set ion mass
  IonMass = ProtonMass();
  
  // revise crossing angle
  crossAng *= 1e-3; // mrad -> rad
  crossAng = -1*TMath::Abs(crossAng); // take -1*abs(crossAng) to enforce the correct sign

  // TEST SETTINGS - for debugging calculations ///////////////
  // set main frame, used for calculations where there is ambiguity which frame is the correct frame to use
  mainFrame = fHeadOn; // fLab, fHeadOn
  // set method for determining `vecQ` 4-momentum components for certain recon methods (JB,DA,(e)Sigma)
  qComponentsMethod = qQuadratic; // qQuadratic, qHadronic, qElectronic
  /////////////////////////////////////////////////////////////
  
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
  RNG = std::make_unique<TRandomMixMax>(91874); // (TODO: fixed seed?)

  // reset counters
  countPIDsmeared=countPIDtrue=0;
};


// ---------------------------------
// calculates 4-momenta components of q and W (`vecQ` and `vecW`) as well as
// derived invariants `W` and `nu`
// - use the quadratic method
// - requires `Q2`, `y`, `Pxh`, `Pyh`
void Kinematics::GetQWNu_quadratic(){

  double f,px,py,pz,pE;
  switch(mainFrame) {
    case fLab:
      f = y*(vecIonBeam.Dot(vecEleBeam));
      pz = vecIonBeam.Pz();
      py = vecIonBeam.Py();
      px = vecIonBeam.Px();
      pE = vecIonBeam.E();
      break;
    case fHeadOn:
      f = y*(HvecIonBeam.Dot(HvecEleBeam));
      pz = HvecIonBeam.Pz();
      py = HvecIonBeam.Py();
      px = HvecIonBeam.Px();
      pE = HvecIonBeam.E();
      break;
  };
  double hx = Pxh - px;
  double hy = Pyh - py;

  double a = 1.0 - (pE*pE)/(pz*pz);
  double b = (2*pE/(pz*pz))*(px*hx + py*hy + f);
  double c = Q2 - hx*hx - hy*hy - (1/(pz*pz))*pow( (f+px*hx+py*hy) ,2.0);
  double disc = b*b - 4*a*c; // discriminant

  double qz1, qz2, qE1, qE2, qE, qz;
  if(disc>=0 && TMath::Abs(a)>1e-6) {
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

    vecQ.SetPxPyPzE(hx, hy, qz, qE);
    if(mainFrame==fHeadOn) 
      this->TransformBackToLabFrame(vecQ,vecQ); // need to be in lab frame for downstream `CalculateHadronKinematics`
    vecW = vecIonBeam + vecQ;
    W = vecW.M();
    Nu = vecIonBeam.Dot(vecQ)/IonMass;

  } else {
    // this happens a lot more often if mainFrame==fLab
    cerr << "ERROR: in Kinematics::GetQWNu_quadratic, ";
    if(isnan(disc))             cerr << "discriminant is NaN";
    else if(disc<0)             cerr << "negative discriminant";
    else if(TMath::Abs(a)<1e-6) cerr << "zero denominator";
    else                        cerr << "unknown reason";
    cerr << "; skipping event" << endl;
    // cerr << "       p=(" << px << "," << py << "," << pz << "," << pE << ") " << endl;
    // cerr << "       a=" << a << endl;
    // cerr << "       b=" << b << endl;
    // cerr << "       c=" << c << endl;
    // cerr << "       disc=" << disc << endl;
    // cerr << "       hx=" << hx << endl;
    // cerr << "       hy=" << hy << endl;
    // cerr << "       Pxh=" << Pxh << endl;
    // cerr << "       Pyh=" << Pyh << endl;
    // cerr << "       f=" << f << endl;
    // cerr << "       y=" << y << endl;
    // cerr << "       Q2=" << Q2 << endl;
    reconOK = false;
  }
};


// calculates 4-momenta components of q and W (`vecQ` and `vecW`) as well as
// derived invariants `W` and `nu`
// - use the hadronic sum 4-momentum `hadronSumVec`
void Kinematics::GetQWNu_hadronic(){
  switch(mainFrame) {
    case fLab:
      vecQ = hadronSumVec - vecIonBeam;
      vecW = hadronSumVec;
      break;
    case fHeadOn:
      vecQ = hadronSumVec - HvecIonBeam;
      vecW = hadronSumVec;
      this->TransformBackToLabFrame(vecQ,vecQ); // need to be in lab frame for downstream `CalculateHadronKinematics`
      this->TransformBackToLabFrame(vecW,vecW);
      break;
  }
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

void Kinematics::GetQWNu_ML(){
  hfsinfo.clear();
  float pidadj = 0;
  if(nHFS >= 2){
    std::vector<float> partHold;
    for(int i = 0; i < nHFS; i++){
      double pidsgn=(hfspid[i]/abs(hfspid[i]));
      if(abs(hfspid[i])==211) pidadj = 0.4*pidsgn;
      if(abs(hfspid[i])==22) pidadj = 0.2*pidsgn;
      if(abs(hfspid[i])==321) pidadj = 0.6*pidsgn;
      if(abs(hfspid[i])==2212) pidadj = 0.8*pidsgn;
      if(abs(hfspid[i])==11) pidadj = 1.0*pidsgn;      
      partHold.push_back(hfspx[i]);
      partHold.push_back(hfspy[i]);
      partHold.push_back(hfspz[i]);
      partHold.push_back(hfsE[i]);
      partHold.push_back(pidadj);
      hfsinfo.push_back(partHold);
      partHold.clear();
    }
    double Q2ele, Q2DA, Q2JB;
    double xele, xDA, xJB;
    TLorentzVector vecQEle;
    globalinfo.clear();
    this->CalculateDISbyElectron();
    vecQEle.SetPxPyPzE(vecQ.Px(), vecQ.Py(), vecQ.Pz(), vecQ.E());
    Q2ele = Q2;
    xele = x;
    this->CalculateDISbyDA();
    Q2DA = Q2;
    xDA = x;
    this->CalculateDISbyJB();
    Q2JB = Q2;
    xJB = x;
    if( Q2DA > 0 && Q2DA < 1e4){
      globalinfo.push_back(log10(Q2DA));
    }
    else{
      globalinfo.push_back(log10((float) (rand()) / (float) (RAND_MAX/10000.0)));
    }
    if( Q2ele > 0 && Q2ele < 1e4){
      globalinfo.push_back(log10(Q2ele));
    }
    else{
      globalinfo.push_back(log10((float) (rand()) / (float) (RAND_MAX/10000.0)));
    }
    if( Q2JB > 0 && Q2JB < 1e4){
      globalinfo.push_back(log10(Q2JB));
    }
    else{
      globalinfo.push_back(log10((float) (rand()) / (float) (RAND_MAX/10000.0)));
    }        
    if(xDA>0 && xDA < 1){
      globalinfo.push_back(-1*log10(xDA));
    }
    else{
      globalinfo.push_back(-1*log10( (float) (rand()) / (float) (RAND_MAX/1.0)  ));
    }
    if(xele>0 && xele < 1){
      globalinfo.push_back(-1*log10(xele));
    }
    else{
      globalinfo.push_back(-1*log10( (float) (rand()) / (float) (RAND_MAX/1.0)  ));
    }
    if(xJB>0 && xJB < 1){
      globalinfo.push_back(-1*log10(xJB));
    }
    else{
      globalinfo.push_back( -1*log10((float) (rand()) / (float) (RAND_MAX/1.0) ) );
    }
    globalinfo.push_back(vecQEle.Px());
    globalinfo.push_back(vecQEle.Py());
    globalinfo.push_back(vecQEle.Pz());
    globalinfo.push_back(vecQEle.E());
#ifdef SIDIS_MLPRED
    py::object nnoutput = pfnimport(hfsinfo, globalinfo);
    std::vector<float> nnvecq = nnoutput.cast<std::vector<float>>();
    vecQ.SetPxPyPzE(nnvecq[0],nnvecq[1],nnvecq[2],nnvecq[3]);
#endif
  }
  else{
    this->CalculateDISbyElectron();
  }
  vecW = vecEleBeam + vecIonBeam - vecElectron; 
  W = vecW.M();
  Nu = vecIonBeam.Dot(vecQ) / IonMass;
  
}

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
  if     (recmethod.CompareTo( "Ele", TString::kIgnoreCase)==0)    { this->CalculateDISbyElectron(); }
  else if(recmethod.CompareTo( "DA", TString::kIgnoreCase)==0)     { this->CalculateDISbyDA(); }
  else if(recmethod.CompareTo( "JB", TString::kIgnoreCase)==0)     { this->CalculateDISbyJB(); }
  else if(recmethod.CompareTo( "Mixed", TString::kIgnoreCase)==0)  { this->CalculateDISbyMixed(); }
  else if(recmethod.CompareTo( "Sigma", TString::kIgnoreCase)==0)  { this->CalculateDISbySigma(); }
  else if(recmethod.CompareTo( "eSigma", TString::kIgnoreCase)==0) { this->CalculateDISbyeSigma(); }
  else if(recmethod.CompareTo( "ML", TString::kIgnoreCase)==0) { this->CalculateDISbyML(); }
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
void Kinematics::CalculateDISbyJB(){
  switch(mainFrame) {
    case fLab:    y = sigmah/(2*vecEleBeam.E()); break;
    case fHeadOn: y = sigmah/(2*HvecEleBeam.E()); break;
  };
  Q2 = (Pxh*Pxh + Pyh*Pyh)/(1-y);
  x = Q2/(s*y);
  switch(qComponentsMethod) {
    case qQuadratic:  this->GetQWNu_quadratic(); break;
    case qHadronic:   this->GetQWNu_hadronic(); break;
    case qElectronic: this->GetQWNu_electronic(); break;
  }
};

// calculate DIS kinematics using DA method
// requires 'vecElectron' set
void Kinematics::CalculateDISbyDA(){
  float thetah = acos( (Pxh*Pxh+Pyh*Pyh - sigmah*sigmah)/(Pxh*Pxh+Pyh*Pyh+sigmah*sigmah) );
  float thetae;
  switch(mainFrame) {
    case fLab:
      thetae = vecElectron.Theta();
      Q2 = 4.0*vecEleBeam.E()*vecEleBeam.E()*sin(thetah)*(1+cos(thetae))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
      break;
    case fHeadOn:
      thetae = HvecElectron.Theta();
      Q2 = 4.0*HvecEleBeam.E()*HvecEleBeam.E()*sin(thetah)*(1+cos(thetae))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
      break;
  };
  y = (sin(thetae)*(1-cos(thetah)))/(sin(thetah)+sin(thetae)-sin(thetah+thetae));
  x = Q2/(s*y);
  switch(qComponentsMethod) {
    case qQuadratic:  this->GetQWNu_quadratic(); break;
    case qHadronic:   this->GetQWNu_hadronic(); break;
    case qElectronic: this->GetQWNu_electronic(); break;
  }
};

// calculate DIS kinematics using mixed method
// requires 'vecElectron' set
void Kinematics::CalculateDISbyMixed(){
  this->GetQWNu_electronic();
  Q2 = -1*vecQ.M2();
  switch(mainFrame) {
    case fLab:    y = sigmah/(2*vecEleBeam.E()); break;
    case fHeadOn: y = sigmah/(2*HvecEleBeam.E()); break;
  }
  x = Q2/(s*y);
};

// calculate DIS kinematics using Sigma method
// requires 'vecElectron' set
void Kinematics::CalculateDISbySigma(){
  switch(mainFrame) {
    case fLab:
      y = sigmah/(sigmah + vecElectron.E()*(1-cos(vecElectron.Theta())));
      Q2 = (vecElectron.Px()*vecElectron.Px() + vecElectron.Py()*vecElectron.Py())/(1-y);
      break;
    case fHeadOn:
      y = sigmah/(sigmah + HvecElectron.E()*(1-cos(HvecElectron.Theta())));
      Q2 = (HvecElectron.Px()*HvecElectron.Px() + HvecElectron.Py()*HvecElectron.Py())/(1-y);
      break;
  }
  x = Q2/(s*y);
  switch(qComponentsMethod) {
    case qQuadratic:  this->GetQWNu_quadratic(); break;
    case qHadronic:   this->GetQWNu_hadronic(); break;
    case qElectronic: this->GetQWNu_electronic(); break;
  }
};

// calculate DIS kinematics using eSigma method
// requires 'vecElectron' set
void Kinematics::CalculateDISbyeSigma(){
  this->CalculateDISbySigma();
  Double_t xsigma = x;
  Double_t ysigma = y; 
  vecQ = vecEleBeam - vecElectron;
  Q2 = -1 * vecQ.M2();
  y = Q2/(s*xsigma);
  x = xsigma;
  switch(qComponentsMethod) {
    case qQuadratic:  this->GetQWNu_quadratic(); break;
    case qHadronic:   this->GetQWNu_hadronic(); break;
    case qElectronic: this->GetQWNu_electronic(); break;
  }
};

void Kinematics::CalculateDISbyML() {
  this->GetQWNu_ML(); // set `vecQ`, `vecW`, `W`, `Nu`   
  Q2 = -1 * vecQ.M2();
  x = Q2 / ( 2 * vecQ.Dot(vecIonBeam) );
  y = vecIonBeam.Dot(vecQ) / vecIonBeam.Dot(vecEleBeam);
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

void Kinematics::CalculateDihadronKinematics() {
  dihadron_phiH = 0;
  dihadron_phiRperp = 1;
  dihadron_theta = 2;
  dihadron_Mh = 3;
  dihadron_pt = 4;
  dihadron_ptLab = 5;
  dihadron_z = 6;
}

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


// add a 4-momentum to the hadronic final state
void Kinematics::AddToHFS(TLorentzVector p4_) {
  TLorentzVector p4 = p4_;
  if(mainFrame==fHeadOn) this->TransformToHeadOnFrame(p4,p4);
  sigmah += (p4.E() - p4.Pz());
  Pxh += p4.Px();
  Pyh += p4.Py();
  hadronSumVec += p4;
  countHadrons++;
};

// add  information from a reconstructed particle to HFSTree
// branches containing full HFS
void Kinematics::AddToHFSTree(TLorentzVector p4, int pid) {
  hfspx.push_back(p4.Px());
  hfspy.push_back(p4.Py());
  hfspz.push_back(p4.Pz());
  hfsE.push_back(p4.E());
  hfspid.push_back(pid);

  nHFS++;
};

// add a specific track and its matched true four momentum
// to HFSTree for final kinematic calculations (not full HFS)
void Kinematics::AddTrackToHFSTree(TLorentzVector p4, int pid) {
  selectedHadPx.push_back(p4.Px());
  selectedHadPy.push_back(p4.Py());
  selectedHadPz.push_back(p4.Pz());
  selectedHadE.push_back(p4.E());
  selectedHadPID.push_back(pid);  
};

// subtract electron from hadronic final state variables
void Kinematics::SubtractElectronFromHFS() {
  if(!isnan(vecElectron.E())){
    switch(mainFrame) {
      case fLab:
        sigmah -= (vecElectron.E() - vecElectron.Pz());
        Pxh -= vecElectron.Px();
        Pyh -= vecElectron.Py();
        hadronSumVec -= vecElectron;
        break;
      case fHeadOn:
        this->TransformToHeadOnFrame(vecElectron,HvecElectron);
        sigmah -= (HvecElectron.E() - HvecElectron.Pz());
        Pxh -= HvecElectron.Px();
        Pyh -= HvecElectron.Py();
        hadronSumVec -= HvecElectron;
        break;
    }
    countHadrons--;    
  } else {
    cerr << "ERROR: electron energy is NaN" << endl;
    // TODO: kill event
  }
};


// reset some variables for the hadronic final state
void Kinematics::ResetHFS() {  
  sigmah = Pxh = Pyh = 0;
  hadronSumVec.SetPxPyPzE(0,0,0,0);
  countHadrons = 0;
  nHFS = 0;
  
  hfspx.clear();
  hfspy.clear();
  hfspz.clear();
  hfsE.clear();
  hfspid.clear();
  
  selectedHadPx.clear();
  selectedHadPy.clear();
  selectedHadPz.clear();
  selectedHadE.clear();
  selectedHadPID.clear();    
};


// DELPHES-only methods //////////////////////////////////////////////////////
#ifndef EXCLUDE_DELPHES
// calculates reconstructed hadronic final state variables from DELPHES tree branches
// expects 'vecElectron' set
// - calculates `sigmah`, `Pxh`, and `Pyh` in the lab and head-on frames
void Kinematics::GetHFS(
    TObjArrayIter itTrack,
    TObjArrayIter itEFlowTrack,
    TObjArrayIter itEFlowPhoton,
    TObjArrayIter itEFlowNeutralHadron,
    TObjArrayIter itpfRICHTrack,
    TObjArrayIter itDIRCepidTrack,   TObjArrayIter itDIRChpidTrack,
    TObjArrayIter itBTOFepidTrack,   TObjArrayIter itBTOFhpidTrack,
    TObjArrayIter itdualRICHagTrack, TObjArrayIter itdualRICHcfTrack
    ) {

  // resets
  this->ResetHFS();
  itTrack.Reset();
  itEFlowTrack.Reset();
  itEFlowPhoton.Reset();
  itEFlowNeutralHadron.Reset();

  // track loop
  while(Track *track = (Track*)itTrack() ){
    TLorentzVector  trackp4 = track->P4();
    if(!isnan(trackp4.E())){
      if( std::abs(track->Eta) < 4.0  ){

        int pid = GetTrackPID( // get smeared PID
            track,
            itpfRICHTrack,
            itDIRCepidTrack, itDIRChpidTrack,
            itBTOFepidTrack, itBTOFhpidTrack,
            itdualRICHagTrack, itdualRICHcfTrack
            );

        if(pid != -1){ // if smeared PID determined, set mass of `trackp4` accordingly
          trackp4.SetPtEtaPhiM(trackp4.Pt(),trackp4.Eta(),trackp4.Phi(),correctMass(pid));
          countPIDsmeared++;
        }
        else { // otherwise, assume true PID
          countPIDtrue++;
          //continue; // drop events if PID not smeared // TODO [critical]: this is more realistic, but resolutions are much worse
        }

        this->AddToHFS(trackp4);
	this->AddToHFSTree(trackp4,pid);
      }
    }    
  }

  // eflow high |eta| track loop
  while(Track *eflowTrack = (Track*)itEFlowTrack() ){
    TLorentzVector eflowTrackp4 = eflowTrack->P4();
    int pid = GetTrackPID( // get smeared PID                                                                     
            eflowTrack,
            itpfRICHTrack,
            itDIRCepidTrack, itDIRChpidTrack,
            itBTOFepidTrack, itBTOFhpidTrack,
            itdualRICHagTrack, itdualRICHcfTrack
            );

    if(!isnan(eflowTrackp4.E())){
      if(std::abs(eflowTrack->Eta) >= 4.0){
        this->AddToHFS(eflowTrackp4);
	this->AddToHFSTree(eflowTrackp4, pid);
      }
    }
  }
  
  // eflow photon loop
  while(Tower* towerPhoton = (Tower*)itEFlowPhoton() ){
    TLorentzVector  towerPhotonp4 = towerPhoton->P4();
    if(!isnan(towerPhotonp4.E())){
      if( std::abs(towerPhoton->Eta) < 4.0  ){
        this->AddToHFS(towerPhotonp4);
	this->AddToHFSTree(towerPhotonp4,22);
      }
    }
  }

  // eflow neutral hadron loop
  while(Tower* towerNeutralHadron = (Tower*)itEFlowNeutralHadron() ){
    TLorentzVector  towerNeutralHadronp4 = towerNeutralHadron->P4();
    if(!isnan(towerNeutralHadronp4.E())){
      if( std::abs(towerNeutralHadron->Eta) < 4.0 ){
        this->AddToHFS(towerNeutralHadronp4);
	this->AddToHFSTree(towerNeutralHadronp4,-1);
      }
    }
  }
  
  // remove electron from hadronic final state
  this->SubtractElectronFromHFS();
};


// calculates generated truth hadronic final state variables from DELPHES tree branches
void Kinematics::GetTrueHFS(TObjArrayIter itParticle){

  // resets
  this->ResetHFS();
  itParticle.Reset();

  // truth loop
  while(GenParticle *partTrue = (GenParticle*)itParticle() ) {
    if(partTrue->Status == 1) this->AddToHFS(partTrue->P4());
  }

  // remove electron from hadronic final state
  this->SubtractElectronFromHFS();
};


// get PID information from PID systems tracks
int Kinematics::GetTrackPID(
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

#endif // ifndef EXCLUDE_DELPHES
// end DELPHES-only methods //////////////////////////////////////////////////////


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

