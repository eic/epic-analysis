R__LOAD_LIBRARY(Largex)

void test_DAG() {
  DAG *D = new DAG();

  D->Print("initial DAG");

  BinSet *BX = new BinSet("X","X"); BX->BuildBins(2,0,2);
  BinSet *BY = new BinSet("Y","Y"); BY->BuildBins(2,0,2);
  BinSet *BZ = new BinSet("Z","Z"); BZ->BuildBins(2,0,2);

  D->AddLayer(BX);
  D->AddLayer(BY);
  D->AddLayer(BZ);

  D->PrintBreadth("full DAG");
  
  D->PrintDepth("depth traversal:");
  D->PrintLeafPaths();
  
  //D->TraverseDepth( D->GetRootNode(), printOp);

};
