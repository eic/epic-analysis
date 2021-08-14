R__LOAD_LIBRARY(Largex)
#include "DAGnode.h"

void test_DAG() {
  DAG *D = new DAG();

  D->AddNode(new DAGnode(tTop,"top0"));

  D->AddEdge("top0",new DAGnode(tBin,"bin0"));
  D->AddEdge("top0",new DAGnode(tBin,"bin1"));
  D->AddEdge("top0",new DAGnode(tBin,"bin2"));

  D->AddNode(new DAGnode(tBottom,"bottom0"));
  D->AddEdge("bin0","bottom0");
  D->AddEdge(D->GetNode("bin1"),"bottom0");
  D->AddEdge("bin2",D->GetNode("bottom0"));

  D->Print();
};
