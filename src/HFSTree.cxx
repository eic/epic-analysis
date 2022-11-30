// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Connor Pecar

#include "HFSTree.h"

ClassImp(HFSTree)

// constructor
HFSTree::HFSTree(TString treeName_, Kinematics *K_, Kinematics *Ktrue_) 
  : treeName(treeName_)
  , K(K_)
  , Ktrue(Ktrue_)
{
  // branch names are set to match `brufit` implementation
  // (see `https://github.com/c-dilks/dispin/tree/master/src`)
  T = new TTree(treeName,treeName);
  T->Branch("vecElectron", &(K->vecElectron));
  T->Branch("vecElectronTrue", &(Ktrue->vecElectron));
  T->Branch("vecEleBeamTrue", &(Ktrue->vecEleBeam));
  T->Branch("vecIonBeamTrue", &(Ktrue->vecIonBeam));
  T->Branch("nHFS", &(K->nHFS), "nHFS/I");
  T->Branch("hfsp4", &(K->hfsp4));
  //T->Branch("hfspy", &(K->hfspy), "hfspy[nHFS]/F");
  //T->Branch("hfspz", &(K->hfspz), "hfspz[nHFS]/F");
  //T->Branch("hfsE", &(K->hfsE), "hfsE[nHFS]/F");
  T->Branch("hfspid", &(K->hfspid), "hfspid[nHFS]/F");
};

// destructor
HFSTree::~HFSTree() {
};
