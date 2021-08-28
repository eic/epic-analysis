R__LOAD_LIBRARY(Largex)

void test_DAG() {
  DAG *D = new DAG();

  D->Print("initial DAG");

  BinSet *BX = new BinSet("X","X"); BX->BuildBins(2,0,2);
  BinSet *BY = new BinSet("Y","Y"); BY->BuildBins(2,0,2);
  BinSet *BZ = new BinSet("Z","Z"); BZ->BuildBins(2,0,2);

  D->AddLayer(BX);
  D->AddLayer(BZ);
  D->AddLayer(BY);

  D->RepatchToLeaf("Z");

  D->PrintBreadth("full DAG");
  D->PrintDepth("depth traversal:");
  D->PrintLeafPaths();


  D->GetRootNode()->StageInboundOp([](){
   cout << "MAIN HEADER" << endl;
  });

  D->GetNode("Z__control")->StageInboundOp([](){
    cout << "Z LOOP HEADER" << endl;
  });
  D->GetNode("Z__control")->StageOutboundOp([](){
    cout << "Z LOOP FOOTER" << endl;
  });

  D->GetLeafNode()->StageInboundOp([](NodePath P){
    cout << "ALGORITHM on ";
    Node::PrintPath(P);
  });
  D->GetLeafNode()->StageOutboundOp([](){
  });

  D->GetRootNode()->StageOutboundOp([](){
    cout << "MAIN FOOTER" << endl;
  });

  D->ExecuteOps();

};
