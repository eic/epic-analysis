/* SimpleTree
 * - provides a simple tree, for common usage in any Analysis or
 *   Analysis derived class; this tree is designed to be compatible
 *   with BruFit for asymmetry analysis
 *   (see `https://github.com/c-dilks/dispin/tree/master/src`)
 */
#ifndef SimpleTree_
#define SimpleTree_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// largex-eic
#include "Kinematics.h"

// ROOT
#include "TTree.h"

class SimpleTree : public TObject
{
  public:
    SimpleTree(TString treeName_, Kinematics *K_);
    ~SimpleTree();

    TTree *GetTree() { return T; };
    Kinematics *GetKinematics() { return K; };
    void FillTree() { T->Fill(); };
    void WriteTree() { T->Write(); };

  private:
    TTree *T;
    Kinematics *K;
    TString treeName;

  ClassDef(SimpleTree,1);
};

#endif
