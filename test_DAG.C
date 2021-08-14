R__LOAD_LIBRARY(Largex)
#include "DAGnode.h"

void test_DAG() {
  DAG *D = new DAG();

  D->AddNode(new DAGnode(tTop,"top0"));

  D->AddNode(new DAGnode(tBin,"l0_b0"));
  D->AddNode(new DAGnode(tBin,"l0_b1"));
  D->AddNode(new DAGnode(tBin,"l0_b2"));

  D->AddNode(new DAGnode(tBin,"l1_b0"));
  D->AddNode(new DAGnode(tBin,"l1_b1"));

  D->AddNode(new DAGnode(tBottom,"bottom0"));

  D->AddEdge("top0","l0_b0");
  D->AddEdge("top0","l0_b1");
  D->AddEdge("top0","l0_b2");

  D->AddEdge("l0_b0","l1_b0");
  D->AddEdge("l0_b0","l1_b1");

  D->AddEdge("l0_b1","l1_b0");
  D->AddEdge("l0_b1","l1_b1");

  D->AddEdge("l0_b2","l1_b0");
  D->AddEdge("l0_b2","l1_b1");

  D->AddEdge("l1_b0","bottom0");
  D->AddEdge("l1_b1","bottom0");

  D->Print();
};
