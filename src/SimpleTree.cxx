#include "SimpleTree.h"

ClassImp(SimpleTree)

// constructor
SimpleTree::SimpleTree(TString treeName_, Kinematics *K_, Kinematics *Ktrue_, int flag=0) 
  : treeName(treeName_)
  , K(K_)
  , Ktrue(Ktrue_)
{
  // branch names are set to match `brufit` implementation
  // (see `https://github.com/c-dilks/dispin/tree/master/src`)
  T = new TTree(treeName,treeName);
  // Original tree from repo
  if(flag==0){
    T->Branch("QSq",       &(K->Q2),        "QSq/D");
    T->Branch("X",         &(K->x),         "X/D");
    T->Branch("Y",         &(K->y),         "Y/D");
    T->Branch("Z",         &(K->z),         "Z/D");
    T->Branch("W",         &(K->W),         "W/D");
    T->Branch("MX",        &(K->mX),        "MX/D");
    T->Branch("PhPerp",    &(K->pT),        "PhPerp/D");
    T->Branch("PhiH",      &(K->phiH),      "PhiH/D");
    T->Branch("PhiS",      &(K->phiS),      "PhiS/D");
    T->Branch("TruePhiH",  &(Ktrue->phiH),  "TruePhiH/D");
    T->Branch("TruePhiS",  &(Ktrue->phiS),  "TruePhiS/D");
    T->Branch("PolT",      &(K->polT),      "PolT/D");
    T->Branch("PolL",      &(K->polL),      "PolL/D");
    T->Branch("PolB",      &(K->polBeam),   "PolB/D");
    T->Branch("Depol1",    &(K->depolP1),   "Depol1/D");
    T->Branch("Depol2",    &(K->depolP2),   "Depol2/D");
    T->Branch("Depol3",    &(K->depolP3),   "Depol3/D");
    T->Branch("Depol4",    &(K->depolP4),   "Depol4/D");
    T->Branch("HadPID",    &(K->hadPID),    "HadPID/I");
    T->Branch("Spin_idx",  &(K->tSpin),     "Spin_idx/I");
    T->Branch("SpinL_idx", &(K->lSpin),     "SpinL_idx/I");
    T->Branch("Weight",    &(weight),       "Weight/D");
  }
  else if(flag==1) // Modified tree for symbolic regression studies
    {
      T->Branch("trueQ2",       &(Ktrue->Q2),        "trueQ2/D");
      T->Branch("trueX",         &(Ktrue->x),         "trueX/D");
      T->Branch("trueY",         &(Ktrue->y),         "trueY/D");
      T->Branch("trues",         &(K->s),             "trues/D");
      
      T->Branch("Q2_e",       &(K->Q2_e),        "Q2_e/D");
      T->Branch("X_e",         &(K->x_e),         "X_e/D");
      T->Branch("Y_e",         &(K->y_e),         "Y_e/D");

      T->Branch("Q2_JB",       &(K->Q2_JB),        "Q2_JB/D");
      T->Branch("X_JB",         &(K->x_JB),         "X_JB/D");
      T->Branch("Y_JB",         &(K->y_JB),         "Y_JB/D");

      T->Branch("Q2_DA",       &(K->Q2_DA),        "Q2_DA/D");
      T->Branch("X_DA",         &(K->x_DA),         "X_DA/D");
      T->Branch("Y_DA",         &(K->y_DA),         "Y_DA/D");

      T->Branch("Ee_i",         &(K->e_Ei),         "e_Ei/D");
      T->Branch("Ee_f",         &(K->e_Ef),         "e_Ef/D");
      T->Branch("thetae",         &(K->e_th),         "e_th/D");

      T->Branch("deltah",         &(K->sigmah),         "deltah/D");
      T->Branch("pxh",         &(K->Pxh),         "Pxh/D");
      T->Branch("pyh",         &(K->Pyh),         "Pyh/D");
      T->Branch("gammah",       &(K->thetah),       "h_th/D");
    }
};
  
// destructor
SimpleTree::~SimpleTree() {
};

