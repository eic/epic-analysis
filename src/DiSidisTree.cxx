// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

#include "DiSidisTree.h"

ClassImp(DiSidisTree)

// constructor
DiSidisTree::DiSidisTree(TString treeName_, std::shared_ptr<Kinematics> K_, std::shared_ptr<Kinematics> Ktrue_) 
  : treeName(treeName_)
  , K(K_)
  , Ktrue(Ktrue_)
{
  // branch names are set to match `brufit` implementation
  // (see `https://github.com/c-dilks/dispin/tree/master/src`)
  T = new TTree(treeName,treeName);
  T->Branch("Q2",            &(K->Q2),                    "Q2/D");
  T->Branch("X",             &(K->x),                     "X/D");
  T->Branch("Y",             &(K->y),                     "Y/D");
  T->Branch("Z",             &(K->dihadron_z),            "Z/D");
  T->Branch("Z1",            &(K->dihadron_z1),           "Z1/D");
  T->Branch("Z2",            &(K->dihadron_z2),           "Z2/D");
  T->Branch("XF",            &(K->dihadron_xF),           "XF/D");
  T->Branch("XF1",           &(K->dihadron_xF1),          "XF1/D");
  T->Branch("XF2",           &(K->dihadron_xF2),          "XF2/D");
  T->Branch("W",             &(K->W),                     "W/D");
  T->Branch("MX",            &(K->dihadron_mX),           "MX/D");
  T->Branch("Mh",            &(K->dihadron_Mh),           "Mh/D");
  T->Branch("PhPerp",        &(K->dihadron_pT),           "PhPerp/D");
  T->Branch("PhPerpLab",     &(K->dihadron_pTlab),        "PhPerpLab/D");
  T->Branch("PhiH",          &(K->dihadron_phiH),         "PhiH/D");
  T->Branch("PhiS",          &(K->phiS),                  "PhiS/D");
  T->Branch("PhiRT",         &(K->dihadron_phiRT),        "PhiRT/D");
  T->Branch("PhiRperp",      &(K->dihadron_phiRperp),     "PhiRperp/D");
  T->Branch("ThetaCOM",      &(K->dihadron_theta),        "ThetaCOM/D");
  T->Branch("TruePhiH",      &(Ktrue->dihadron_phiH),     "TruePhiH/D");
  T->Branch("TruePhiS",      &(Ktrue->phiS),              "TruePhiS/D");
  T->Branch("TruePhiRT",     &(Ktrue->dihadron_phiRT),    "TruePhiRT/D");
  T->Branch("TruePhiRperp",  &(Ktrue->dihadron_phiRperp), "TruePhiRperp/D");
  T->Branch("TrueThetaCOM",  &(Ktrue->dihadron_theta),    "TrueThetaCOM/D");
  T->Branch("PolT",          &(K->polT),                  "PolT/D");
  T->Branch("PolL",          &(K->polL),                  "PolL/D");
  T->Branch("PolB",          &(K->polBeam),               "PolB/D");
  T->Branch("Spin_idx",      &(K->tSpin),                 "Spin_idx/I");
  T->Branch("SpinL_idx",     &(K->lSpin),                 "SpinL_idx/I");
  T->Branch("Weight",        &(weight),                   "Weight/D");
};

// destructor
DiSidisTree::~DiSidisTree() {
};

