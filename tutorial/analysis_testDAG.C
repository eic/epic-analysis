R__LOAD_LIBRARY(Largex)

// test DAG implementation
void analysis_testDAG(
    TString infiles="datarec/tutorial.config", /* list of input files */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="testDAG" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 100;

  // set reconstruction method ====================================
  A->SetReconMethod("Ele");

  // set binning scheme ====================================

  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBins(2,0.03,1,true);
  
  A->AddBinScheme("x");
  A->BinScheme("x")->BuildBins(2,0.05,1,true);

  A->AddBinScheme("q2");
  A->BinScheme("q2")->BuildBins(2,1,100,true);

  /*
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Max",0.1);
  A->BinScheme("y")->BuildBin("Min",0.1);
  */


  // perform the analysis ==================================
  A->Execute();
};
