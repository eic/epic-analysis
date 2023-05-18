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

  T->Branch("selectedHadPx", &(K->selectedHadPx));
  T->Branch("selectedHadPy", &(K->selectedHadPy));
  T->Branch("selectedHadPz", &(K->selectedHadPz));
  T->Branch("selectedHadE", &(K->selectedHadE));
  T->Branch("selectedHadPID", &(K->selectedHadPID));

  T->Branch("selectedHadPxTrue", &(Ktrue->selectedHadPx));
  T->Branch("selectedHadPyTrue", &(Ktrue->selectedHadPy));
  T->Branch("selectedHadPzTrue", &(Ktrue->selectedHadPz));
  T->Branch("selectedHadETrue", &(Ktrue->selectedHadE));
  T->Branch("selectedHadPIDTrue", &(Ktrue->selectedHadPID));
  
  T->Branch("weight",    &(weight),       "weight/D");
};

// destructor
HFSTree::~HFSTree() {
};

