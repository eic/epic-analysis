// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek

/* DiSidisTree
 * - provides a dihadron tree
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

class DiSidisTree
{
  public:
    DiSidisTree(TString treeName_, std::shared_ptr<Kinematics> K_, std::shared_ptr<Kinematics> Ktrue_);
    ~DiSidisTree();

    TTree *GetTree() { return T; };
    std::shared_ptr<Kinematics> GetKinematics() { return K; };
    std::shared_ptr<Kinematics> GetKinematicsTrue() { return Ktrue; };
    void FillTree(Double_t w) { weight = w; T->Fill(); };
    void WriteTree() { T->Write(); };

  private:
    Double_t weight;
    TTree *T;
    std::shared_ptr<Kinematics> K;
    std::shared_ptr<Kinematics> Ktrue;
    TString treeName;

  ClassDef(DiSidisTree,1);
};
