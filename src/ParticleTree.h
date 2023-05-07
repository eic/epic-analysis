// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

/* ParticleTree
 * - provides a particle tree for storing reconstructed particle kinematics + pid
 * - helpful for debugging matching algorithm
 */
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// epic-analysis

// ROOT
#include <TTree.h>
#include <TLorentzVector.h>

class ParticleTree
{
  public:
    ParticleTree(TString treeName_);
    ~ParticleTree();

    TTree *GetTree() { return T; };
    void FillTree(TLorentzVector recopart, TLorentzVector mcpart, int pid, int status) {
      recopart_ = recopart;
      mcpart_   = mcpart;
      pid_      = pid;
      status_   = status;
      T->Fill(); };
    void WriteTree() { T->Write(); };

  private:
  
    TTree *T;
    TString treeName;
    TLorentzVector recopart_;
    TLorentzVector mcpart_;
    int status_;
    int pid_;
  
  ClassDef(ParticleTree,1);
};
