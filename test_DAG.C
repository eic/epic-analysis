R__LOAD_LIBRARY(Largex)
#include "DAGnode.h"

void test_DAG() {
  DAG *D = new DAG();

  D->Print("initial DAG");

  std::vector<DAGnode*> lay0;
  lay0.push_back(new DAGnode(tBin,"l0_b0"));
  lay0.push_back(new DAGnode(tBin,"l0_b1"));
  lay0.push_back(new DAGnode(tBin,"l0_b2"));

  D->AddLayer(lay0);

  D->Print("final DAG");
};
