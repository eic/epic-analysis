// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make tables of kinematics averages
void postprocess_qbins_tables(
    TString infile="out/tutorial.qbins.root"
) {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // print DAG ==================================================
  P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // start output file  ========================================
  // - we just store it in the "images" directory
  TString tableFile = P->GetPngDir() + "/table_counts_qbins.txt";
  P->StartTextFile(tableFile,"Counts and averages in bins of Q2");


  // lambdas ===================================================
  
  // before Q2 subloop operator:
  // print a header for the current bin; to be run before Q2 bin loop
  // - we must capture `P`, the pointer to `PostProcessor`, along
  //   with the name of our output file
  // - use the `NodePath` to print the bin number; the `NodePath` will
  //   give the known bin numbers just before the Q2 bin loop
  // - see `../src/NodePath.h` for the methods
  // - `NodePath` is a linked list of the current path in the DAG during
  //   depth-first traversal
  auto printHeader = [&P,&tableFile](NodePath *bins){
    // append a line to the specified text file
    P->AppendToTextFile(
        tableFile,
        "\nKinematic Bin: " + bins->GetBinNode("pt")->GetID() + 
        "\n=========================================================="
        );
  };

  // payload operator:
  // - this lambda expression will be executed on every multidimensional
  //   bin's Histos object; in our case it will run on every Q2 bin
  //   (nested in the loop over pT bins)
  // - it will print out a row of average kinematic values for the 
  //   current Q2 bin, and append it to the table file
  auto payloadOp = [&P,&tableFile](Histos *H){
    // PostProcessor::DumpAve dumps a table of average values,
    // in bins of Q2; each row will be for one bin
    P->DumpAve( tableFile, H, "q2");
  };

  // after Q2 subloop operator:
  // - finish algorithm, called after loop over table rows, so that
  //   PostProcessor knows to start a new table for the next iteration
  // - we do not need any lambda arguments
  auto finishTable = [&P,&tableFile](){
    P->FinishDumpAve(tableFile);
  };


  // stage lambdas ===================================================
  // - now we put ("stage") these lambdas on the DAG, by creating a subloop
  // - the first argument is a list of bin layers to include in the 
  //   subloop, in the order you want; a control node will be created
  //   before the first specified layer; in this case we only make
  //   a subloop over Q2 bins
  // - the second argument is the "inbound" lambda; it will be executed
  //   before looping over Q2; this argument is optional
  // - the third argument is the "outbound" lambda; it will be executed
  //   after looping over Q2; this argument is optional
  P->Op()->Subloop(
      {"q2"}, // layers, given as a std::vector<TString>
      printHeader, // inbound lambda (executed before subloop)
      finishTable // outbound lambda (executed after subloop)
      );
  // - then we stage the payload lambda, which runs on every Q2 bin
  P->Op()->Payload(payloadOp);


  // print DAG ============================
  // - compare the configured DAG with the initial DAG; notice
  //   the new control node before the Q2 bins; you have staged the lambdas
  //   `printHeader` and `finishTable` on that node
  // - you have also staged `payloadOp` on the leaf node
  // - the lambdas you have staged are not printed; only the node names
  //   and their connections (inputs/outputs) are printed
  P->Op()->PrintBreadth("HistosDAG Final Setup");

  // execution ===================================================
  // - after you have defined all your operators, run them:
  P->Execute();
  
  // finish ===================================================
  // dump the table to stdout
  P->PrintTextFile(tableFile);
  cout << tableFile << " written" << endl;
  // - PostProcessor::Finish() is NECESSARY to close output streams
  // - output file names and directories will be printed to stdout
  P->Finish(); 
};
