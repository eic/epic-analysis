R__LOAD_LIBRARY(Largex)

// make tables of kinematics averages
void postprocess_qbins_tables(
    TString infile="out/tutorial.xqbins.example_5x41.root"
) {


  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  Analysis *A = new Analysis();

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // start output file  ========================================
  // - we just store it in the "images" directory
  TString tableFile = P->GetPngDir() + "/table_counts_qbins.txt";
  P->StartTextFile(tableFile,"Counts and averages in bins of Q2");


  // loop over bins ==================================================
  for(int bx  : P->GetBinNums("x")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int by  : P->GetBinNums("y")) {
  for(int bfs : P->GetBinNums("finalState")) {

    // loop over pT bins
    for(int bpt : P->GetBinNums("pt")) {

      //if( P->SkipFull("pt",bpt) ) continue; // skip full bin

      // print bin header
      P->AppendToTextFile(
          tableFile,
          Form("\nKinematic Bin: %d =========================",bpt));

      // loop over Q2 bins
      for(int bq : P->GetBinNums("q2")) {

        if( P->SkipFull("q2",bq) ) continue; // skip full bin

        // PostProcessor::DumpAve dumps a table of average values,
        // in bins of Q2; each row will be for one bin
        P->DumpAve(
            tableFile,
            A->GetHistosName(bpt,bx,bz,bq,by,bfs),
            "q2");

      }; // end loop over Q2 bins

      
      // finish algorithm, called after loop over table rows, so that
      // PostProcessor knows to start a new table for the next iteration
      P->FinishDumpAve(tableFile);

    }; // end pT loop
  }}}}; // end outer binning loop


  // finish ===================================================

  // dump the table to stdout
  P->PrintTextFile(tableFile);
  cout << tableFile << " written" << endl;

  // PostProcessor::Finish() is NECESSARY to close output streams
  // - output file names and directories will be printed to stdout
  P->Finish(); 

};
