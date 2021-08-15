R__LOAD_LIBRARY(Largex)

void test_DAG() {
  DAG *D = new DAG();

  D->Print("initial DAG");

  std::vector<Node*> lay0;

  BinSet *BX = new BinSet("X","X"); BX->BuildBins(2,0,2);
  BinSet *BY = new BinSet("Y","Y"); BY->BuildBins(2,0,2);
  BinSet *BZ = new BinSet("Z","Z"); BZ->BuildBins(2,0,2);

  D->AddLayer(BX);
  D->AddLayer(BY);
  D->AddLayer(BZ);

  D->Print("full DAG");

  D->RepatchToLeaf("Y");
  D->RepatchToLeaf("X");

  //D->RepatchAllToFull();
  //D->RepatchToFull("Y__control");

  D->Print("repatched DAG");
};
