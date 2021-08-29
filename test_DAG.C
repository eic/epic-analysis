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

  D->PrintBreadth("full DAG");
  D->PrintDepth("depth traversal:");
  D->PrintLeafPaths();
  
  // -------------

  D->Initial([](){ cout << "MAIN HEADER" << endl; });
  D->Final([](){ cout << "MAIN FOOTER" << endl; });

  D->Control(
    {"Z","Y"},
    [](){ cout << "Z LOOP HEADER - CONTROL" << endl; },
    [](){ cout << "Z LOOP FOOTER - CONTROL" << endl; }
    );

  //D->Before( {"Z"}, [](){ cout << "Z LOOP HEADER - BEFORE" << endl; });
  //D->After( {"Z"}, [](){ cout << "Z LOOP FOOTER - AFTER" << endl; });
  D->Before({"Y"},[](){  cout << "Y LOOP HEADER" << endl; });

  D->ForEach([](Histos*){}); // TODO

  D->ExecuteOps();

};
