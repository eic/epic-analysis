R__LOAD_LIBRARY(Largex)

/* run in a grid of (x,Q2) 2D bins
 * - various ways to make a grid are demonstrated
 * - observe how the resulting histograms differ in each (x,Q2) bin
 */
void analysis_coverage(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="coverage" /* output filename prefix*/
) {

  // setup analysis ========================================
  Analysis *A = new Analysis(
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
  A->AddBinScheme("x");
  A->AddBinScheme("z");

  A->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );
  A->BinScheme("z")->BuildBin("Range", 0.2, 0.5 );
  A->BinScheme("z")->BuildBin("Range", 0.5, 0.8 );


  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
