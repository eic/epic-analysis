// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Connor Pecar

#include "HFSTree.h"

ClassImp(HFSTree)

// constructor
HFSTree::HFSTree(TString treeName_, std::shared_ptr<Kinematics> K_, std::shared_ptr<Kinematics> Ktrue_) 
  : treeName(treeName_)
  , K(K_)
  , Ktrue(Ktrue_)
{
  T = new TTree(treeName,treeName);
  T->Branch("vecElectron", &(K->vecElectron));
  T->Branch("vecElectronTrue", &(Ktrue->vecElectron));
  T->Branch("vecEleBeamTrue", &(Ktrue->vecEleBeam));
  T->Branch("vecIonBeamTrue", &(Ktrue->vecIonBeam));
  T->Branch("nHFS", &(K->nHFS), "nHFS/I");

  T->Branch("hfspx", &(K->hfspx));
  T->Branch("hfspy", &(K->hfspy));
  T->Branch("hfspz", &(K->hfspz));
  T->Branch("hfsE", &(K->hfsE));
  T->Branch("hfspid", &(K->hfspid));

  T->Branch("hfspxTrue", &(K->hfspxTrue));
  T->Branch("hfspyTrue", &(K->hfspyTrue));
  T->Branch("hfspzTrue", &(K->hfspzTrue));
  T->Branch("hfsETrue", &(K->hfsETrue));
  T->Branch("hfspidTrue", &(K->hfspidTrue));
  
  T->Branch("weight",    &(weight),       "weight/D");
};

// destructor
HFSTree::~HFSTree() {
};

