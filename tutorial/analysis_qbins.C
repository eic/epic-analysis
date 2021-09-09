R__LOAD_LIBRARY(Largex)

// run in Q2 bins, for two pT ranges
void analysis_qbins(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="tutorial.qbins" /* output filename prefix*/
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
  A->BinScheme("q2")->BuildBins( 10, 1, 100, true );

  A->AddBinScheme("pt");
  A->BinScheme("pt")->BuildBin( "Max", 0.5 ); // pT<0.5 GeV
  A->BinScheme("pt")->BuildBin( "Min", 0.5 ); // pT>0.5 GeV

  // perform the analysis ==================================
  A->Execute();

  // for reference, here is a print out of HistosDAG
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
