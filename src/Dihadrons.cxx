// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

#include "Dihadrons.h"

ClassImp(DihadronSet)
ClassImp(Dihadron)

Dihadron::Dihadron(TLorentzVector vecHad1, TLorentzVector vecHad2, Kinematics *K)
{
  const TLorentzVector &vecE   = K->vecEleBeam;
  const TLorentzVector &vecP   = K->vecIonBeam;
  const TLorentzVector &vecL   = K->vecElectron;
  const TLorentzVector &vecQ   = K->vecQ;
  const TLorentzVector &vecW   = K->vecW;
  const TLorentzVector &vecPh  = vecHad1 + vecHad2;
  const TLorentzVector &vecR   = 0.5 * (vecHad1-vecHad2);
  const TLorentzVector vecH[2] = {vecHad1, vecHad2};

  // boost to ion rest frame
  TLorentzVector IvecH[2];
  TLorentzVector IvecL, IvecQ, IvecPh;
  K->BoostToIonFrame(vecL,IvecL);
  K->BoostToIonFrame(vecQ,IvecQ);
  K->BoostToIonFrame(vecPh,IvecPh);
  for(int h=0; h<2; h++)
    K->BoostToIonFrame(vecH[h],IvecH[h]);

  // invariant mass, missing mass, Z, PhPerp, Zeta
  Mh = vecPh.M();
  MX = TMath::Abs((vecW-vecPh).M());
  Z  = vecP.Dot(vecPh) / vecP.Dot(vecQ);
  for(int h=0; h<2; h++)
    hadZ[h] = vecP.Dot(vecH[h]) / vecP.Dot(vecQ);
  Ph     = vecPh.Vect().Mag();
  PhPerp = Kinematics::Reject(IvecPh.Vect(),IvecQ.Vect()).Mag();
  Zeta   = 2 * vecR.Dot(vecP) / vecPh.Dot(vecP);

  // PhiH
  PhiH = Kinematics::PlaneAngle(
      IvecQ.Vect(), IvecL.Vect(),
      IvecQ.Vect(), IvecPh.Vect()
      );

  // PhiR
  Double_t coeff =
    K->x * (Zeta*Mh*Mh - (vecH[0].M2() - vecH[1].M2()) ) /
    ( K->Q2 * Z );
  TLorentzVector vecRT = vecR - (Zeta/2)*vecPh + coeff*vecP;
  TLorentzVector IvecRT;
  K->BoostToIonFrame(vecRT,IvecRT);
  PhiR = Kinematics::PlaneAngle(
      IvecQ.Vect(), IvecL.Vect(),
      IvecQ.Vect(), IvecRT.Vect()
      );

  // PhiS
  PhiS = K->phiS;

  // Theta
  Double_t MRterm[2];
  for(int h=0; h<2; h++)
    MRterm[h] = TMath::Sqrt( vecH[h].M2() + vecR.Vect().Mag2() );
  Theta = TMath::ACos(
      ( MRterm[0] - MRterm[1] - Mh*Zeta ) / 
      ( 2 * vecR.Vect().Mag() )  
      );
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

void DihadronSet::IncludeHadron(TString hadName) {
  includedHadrons.push_back(hadName);
}

void DihadronSet::AddHadron(Analysis *A) {
  if(hadSetRec.find(A->finalStateID)==hadSetRec.end()) {
    hadSetRec.insert({A->finalStateID,{}});
    hadSetGen.insert({A->finalStateID,{}});
  }
  if(debug) fmt::print("DihadronSet: AddHadron '{}'\n",A->finalStateID);
  hadSetRec[A->finalStateID].push_back( A->kin->vecHadron     );
  hadSetGen[A->finalStateID].push_back( A->kinTrue->vecHadron );
}

void DihadronSet::CalculateKinematics(Analysis *A, Double_t wgt) {

  // TODO: this does not yet handle generalized dihadrons
  if(includedHadrons.size()!=2) fmt::print("ERROR: more or less than 2 final states defined for DihadronSet\n");

  // hadron pairing
  auto PairHadrons = [this] (auto hadSet, auto &dihList, Kinematics *K) {
    TString hadNames[2] = { includedHadrons[0], includedHadrons[1] };
    for(const auto &had0vec : hadSet[hadNames[0]]) {
      for(const auto &had1vec : hadSet[hadNames[1]]) {
        dihList.push_back(Dihadron(had0vec, had1vec, K));
        if(debug) fmt::print(" pair: {}, {}\n",hadNames[0],hadNames[1]);
      }
    }
  };
  if(debug) fmt::print("DihadronSet::CalculateKinematics reconstructed\n");
  PairHadrons(hadSetRec,dihListRec,A->kin);
  if(debug) fmt::print("DihadronSet::CalculateKinematics generated\n");
  PairHadrons(hadSetGen,dihListGen,A->kinTrue);

  // fill output data structures
  for(std::size_t i=0; i<dihListRec.size(); i++) {
    A->dih     = &(dihListRec[i]);
    A->dihTrue = &(dihListGen[i]);
    auto finalStateID_tmp = A->finalStateID; // temporarily change Analysis::finalStateID
    A->finalStateID = dihadronFinalStateID;  // to that of this DihadronSet
    A->FillHistos2h(wgt);
    A->finalStateID = finalStateID_tmp; // revert Analysis::finalStateID
  }

  // reset internal storage
  hadSetRec.clear();
  hadSetGen.clear();
  dihListRec.clear();
  dihListGen.clear();
}
