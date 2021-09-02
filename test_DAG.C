R__LOAD_LIBRARY(Largex)

void test_DAG() {

  Analysis *A = new Analysis();
  A->AddBinScheme("x"); A->BinScheme("x")->BuildBins(2,0,2);
  A->AddBinScheme("y"); A->BinScheme("y")->BuildBins(2,0,2);
  A->AddBinScheme("z"); A->BinScheme("z")->BuildBins(2,0,2);
  HistosDAG *D = new HistosDAG();
  D->Build(A->GetBinSchemes());


  D->PrintBreadth("breadth traversal:");
  D->PrintDepth("depth traversal:");
  D->PrintLeafPaths();
  
  // -------------

  D->Initial([](){ cout << "MAIN HEADER" << endl; });
  D->Final([](){ cout << "MAIN FOOTER" << endl; });

  /*
  D->Subloop(
    {"Z","Y"},
    [](){ cout << "Z LOOP HEADER - CONTROL" << endl; },
    [](){ cout << "Z LOOP FOOTER - CONTROL" << endl; }
    );
  */

  //D->BeforeSubloop( {"Z"}, [](){ cout << "Z LOOP HEADER - BEFORE" << endl; });
  //D->AfterSubloop( {"Z"}, [](){ cout << "Z LOOP FOOTER - AFTER" << endl; });
  //D->BeforeSubloop({"Y"},[](){  cout << "Y LOOP HEADER" << endl; });

  D->ForEach([](Histos* H){
    std::cout << "PAYLOAD: " << std::endl
              << "   " << H->GetSetName() << std::endl
              << "   " << H->GetSetTitle() << std::endl;
  });

  D->ExecuteOps();

};
