R__LOAD_LIBRARY(Largex)

// make resolution plots
void analysis_resolution(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="resolution" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state
  //A->AddFinalState("jet"); // jets

  // set binning scheme ====================================

  A->AddBinScheme("q2");
  A->BinScheme("q2")->BuildBins( 4, 1,    100,  true );

  A->AddBinScheme("x");
  A->BinScheme("x")->BuildBins(  6, 0.01, 1,    true );

  /*
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
  A->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );
  */


  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
