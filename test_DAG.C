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

  D->GetNode("Z__control")->StageInboundOp(
    [](Node *N,NodePath P){
      cout << "Z LOOP HEADER" << endl;
    });
  D->GetNode("Z__control")->StageOutboundOp(
    [](Node *N,NodePath P){
      cout << "Z LOOP FOOTER" << endl;
    });

  D->GetLeafNode()->StageInboundOp(
    [&D](Node *N,NodePath P){
      cout << "ALGORITHM on ";
      Node::PrintPath(P);
    });
  D->GetLeafNode()->StageOutboundOp([](Node *N,NodePath P){});

  D->GetRootNode()->StageInboundOp(
    [](Node *N,NodePath P){ cout << "MAIN HEADER" << endl; });
  D->GetRootNode()->StageOutboundOp(
    [](Node *N,NodePath P){ cout << "MAIN FOOTER" << endl; });

  D->ExecuteOps();

};
