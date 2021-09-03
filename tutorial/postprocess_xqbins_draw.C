R__LOAD_LIBRARY(Largex)

// make kinematics coverage plots, such as eta vs. p in bins of (x,Q2)
void postprocess_xqbins_draw(
    TString infile="out/tutorial.xqbins.example_5x41.root"
) {


  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  Analysis *A = new Analysis();

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);


  // loop over bins ==================================================
  // TODO: when HistosDAG is implemented, loops like this will
  // be significantly simpler; for now we unfortunately have to use
  // these large, nested for loops
  for(int bpt : P->GetBinNums("pt")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int by  : P->GetBinNums("y")) {
  for(int bfs : P->GetBinNums("finalState")) {


    /* if you have other bins, such as z-bins, you could print
     * out some information here about which z-bin you are in
     */


    // loop over (x,Q2) grid
    for(int bx  : P->GetBinNums("x")) {
    for(int bq  : P->GetBinNums("q2")) {


      // skip the full-range bins
      // TODO: with HistosDAG this won't be necessary anymore, since
      // full-range bins won't be needed by default (you can still
      // add them if you want)
      if( P->SkipFull("x",bx) ) continue;
      if( P->SkipFull("q2",bq) ) continue;
      if( P->SkipFull("z",bz) ) continue;


      /* PostProcessor::DrawSingle: draw a single plot for this (x,Q2) bin
       * - png files will be generated 
       * - the TCanvas will also be written to a root file
       */

      // Q2 vs. x (sanity check)
      P->DrawSingle(
          A->GetHistosName(bpt,bx,bz,bq,by,bfs),
          "Q2vsX"
          );

      // eta vs. p
      P->DrawSingle(
          A->GetHistosName(bpt,bx,bz,bq,by,bfs),
          "etaVsP"
          );

      /* - see `../src/Analysis.cxx` or the ROOT file for other histograms
       *   that you can draw; you are welcome to add your own
       * - see also `../src/PostProcessor.cxx` for other post-processing
       *   methods; you are welcome to add your own
       */
      

    }}; // end (x,Q2) loop

  }}}}; // end outer binning loop


  // finish ===================================================

  // PostProcessor::Finish() is NECESSARY to close output streams
  // - output file names and directories will be printed to stdout
  P->Finish(); 

};
