// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Duane Byer

/* SidisTree
 * - provides a simple tree, for common usage in any Analysis or
 *   Analysis derived class; this tree is designed to be compatible
 *   with BruFit for asymmetry analysis
 *   (see `https://github.com/c-dilks/dispin/tree/master/src`)
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

class SidisTree
{
  public:
    SidisTree(TString treeName_, std::shared_ptr<Kinematics> K_, std::shared_ptr<Kinematics> Ktrue_);
    ~SidisTree();

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

  ClassDef(SidisTree,1);
};
