// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Connor Pecar

/* SimpleTree
 * - provides a simple tree, for common usage in any Analysis or
 *   Analysis derived class; this tree is designed to be compatible
 *   with BruFit for asymmetry analysis
 *   (see `https://github.com/c-dilks/dispin/tree/master/src`)
 */
#ifndef HFSTree_
#define HFSTree_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// sidis-eic
#include "Kinematics.h"

// ROOT
#include "TTree.h"

class HFSTree : public TObject
{
  public:
    HFSTree(TString treeName_, Kinematics *K_, Kinematics *Ktrue_);
    ~HFSTree();

    TTree *GetTree() { return T; };
    Kinematics *GetKinematics() { return K; };
    Kinematics *GetKinematicsTrue() { return Ktrue; };
    void FillTree(Double_t w) { weight = w;
      T->Fill(); };
    void WriteTree() { T->Write(); };
  
  private:
    Double_t weight;
    TTree *T;
    Kinematics *K;
    Kinematics *Ktrue;
    TString treeName;

  ClassDef(HFSTree,1);
};

#endif
