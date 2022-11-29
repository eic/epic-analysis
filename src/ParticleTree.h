// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Gregory Matousek

/* ParticleTree
 * - provides a particle tree for storing reconstructed particle kinematics + pid
 * - helpful for debugging matching algorithm
 */
#ifndef ParticleTree_
#define ParticleTree_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// sidis-eic

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

class ParticleTree : public TObject
{
  public:
    ParticleTree(TString treeName_);
    ~ParticleTree();

    TTree *GetTree() { return T; };
    void FillTree(TLorentzVector recopart, TLorentzVector mcpart, int pid, int status, Double_t w) {
      recopart_ = recopart;
      mcpart_   = mcpart;
      pid_      = pid;
      status_   = status;
      weight    = w;
      T->Fill(); };
    void WriteTree() { T->Write(); };

  private:
  
    Double_t weight;
    TTree *T;
    TString treeName;
    TLorentzVector recopart_;
    TLorentzVector mcpart_;
    int status_;
    int pid_;
  
  ClassDef(ParticleTree,1);
};

#endif
