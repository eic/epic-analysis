#include "SimpleTree.h"

ClassImp(SimpleTree)

// constructor
SimpleTree::SimpleTree(TString treeName_, Kinematics *K_) 
  : treeName(treeName_)
  , K(K_)
{
  T = new TTree(treeName,treeName);
  T->Branch("x",&(K->x),"x/D");
  T->Branch("Q2",&(K->Q2),"Q2/D");
};

// destructor
SimpleTree::~SimpleTree() {
};

