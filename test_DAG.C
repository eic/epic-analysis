R__LOAD_LIBRARY(Largex)

void test_DAG() {
  DAG *D = new DAG();

  D->Print("initial DAG");

  std::vector<Node*> lay0;
  lay0.push_back(new Node(NT::bin,"l0_b0"));
  lay0.push_back(new Node(NT::bin,"l0_b1"));
  lay0.push_back(new Node(NT::bin,"l0_b2"));

  D->AddLayer(lay0);

  D->Print("final DAG");
};
