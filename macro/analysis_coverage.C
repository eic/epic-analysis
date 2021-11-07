R__LOAD_LIBRARY(Largex)

/* run in a grid of (x,Q2) 2D bins
 * - various ways to make a grid are demonstrated
 * - observe how the resulting histograms differ in each (x,Q2) bin
 */
void analysis_coverage(
    TString infiles="datarec/in.config", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=25, /* crossing angle [mrad] */
    TString outfilePrefix="coverage" /* output filename prefix*/
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
  A->AddFinalState("KpTrack"); // kaon final state
  A->AddFinalState("pTrack"); // proton final state                                                                                                                         
  //A->AddFinalState("jet"); // jets


  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================

  /* TODO
   * - finer binning
   * - different sqrt(s) values
   * - eta vs. p in bins of (x,Q2)
   * - Q2 vs. x in bins of (eta,p)
   */

  A->AddBinScheme("q2");
  A->AddBinScheme("x");
  A->AddBinScheme("z");

  A->BinScheme("q2")->BuildBins( 5, 1,    1000,  true );
  A->BinScheme("x")->BuildBins(  10, 1e-4, 1,    true );
  A->BinScheme("z")->BuildBin("Range", 0.2, 0.4 );
  A->BinScheme("z")->BuildBin("Range", 0.4, 0.8 );
  A->BinScheme("z")->BuildBin("Min", 0.2);


  // perform the analysis ==================================
  A->Execute();

  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
