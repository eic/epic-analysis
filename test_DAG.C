R__LOAD_LIBRARY(Largex)

void test_DAG() {
  HistosDAG *D = new HistosDAG();

  D->Print("initial DAG");

  BinSet *BX = new BinSet("X","X"); BX->BuildBins(2,0,2);
  BinSet *BY = new BinSet("Y","Y"); BY->BuildBins(2,0,2);
  BinSet *BZ = new BinSet("Z","Z"); BZ->BuildBins(2,0,2);

  D->AddLayer(BX);
  D->AddLayer(BZ);
  D->AddLayer(BY);

  D->Control(
    {"Z"},
    [](){ cout << "Z LOOP HEADER" << endl; },
    [](){ cout << "Z LOOP FOOTER" << endl; }
    );

  D->PrintBreadth("full DAG");
  D->PrintDepth("depth traversal:");
  D->PrintLeafPaths();

  D->Initial([](){ cout << "MAIN HEADER" << endl; });
  D->Final([](){ cout << "MAIN FOOTER" << endl; });

  D->Payload([](Histos*){}); // TODO

  D->Execute();

};
