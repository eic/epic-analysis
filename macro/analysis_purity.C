R__LOAD_LIBRARY(Largex)

// ratios of histograms with y-cut enabled to those with y-cut disabled
void analysis_purity(
    TString infiles="datarec/example_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString methodname="ele", /*reconstruction method name*/
    TString outfilePrefix="resolutions" /* output filename prefix*/
    
) {

  //outfilePrefix+="_DA";
  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix+"_"+methodname
      );

  //  A->maxEvents = 100; // use this to limit the number of events
  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  A->SetReconMethod(methodname); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state
  //A->AddFinalState("jet"); // jets

  // set binning scheme ====================================
  // z ranges
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBin("Min",0.2); // needed?

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Full"); // a bin with no y-cut

  // perform the analysis ==================================
  A->Execute();

};
