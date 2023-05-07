// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Connor Pecar

/* HFSTree
   - Produces a tree containing information needed for kinematic
     reconstruction studies: hadronic final state four-momenta
     and PID, scattered electron information, and true scattered
     electron and beam information.
 */
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// epic-analysis
#include "Kinematics.h"

// ROOT
#include <TTree.h>

class HFSTree
{
  public:
    HFSTree(TString treeName_, std::shared_ptr<Kinematics> K_, std::shared_ptr<Kinematics> Ktrue_);
    ~HFSTree();

    TTree *GetTree() { return T; };
    std::shared_ptr<Kinematics> GetKinematics() { return K; };
    std::shared_ptr<Kinematics> GetKinematicsTrue() { return Ktrue; };
    void FillTree(Double_t w) { weight = w;
      T->Fill(); };
    void WriteTree() { T->Write(); };
  
  private:
    Double_t weight;
    TTree *T;
    std::shared_ptr<Kinematics> K;
    std::shared_ptr<Kinematics> Ktrue;
    TString treeName;

  ClassDef(HFSTree,1);
};
