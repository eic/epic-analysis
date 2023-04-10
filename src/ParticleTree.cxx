// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

#include "ParticleTree.h"

ClassImp(ParticleTree)

// constructor
ParticleTree::ParticleTree(TString treeName_) 
  : treeName(treeName_)
{
  T = new TTree(treeName,treeName);
  T->Branch("recPart",  "TLorentzVector" , &(recopart_));    
  T->Branch("mcPart",   "TLorentzVector" , &(mcpart_));    
  T->Branch("pid",      &(pid_)           , "pid/I");
  T->Branch("status",      &(status_)           , "status/I");
};


// destructor
ParticleTree::~ParticleTree() {
};

