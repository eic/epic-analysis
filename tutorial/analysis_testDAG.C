R__LOAD_LIBRARY(Largex)

// test DAG implementation
void analysis_testDAG(
    TString infiles="datarec/tutorial.config", /* list of input files */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=-25, /* crossing angle [mrad] */
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


  // define cuts ====================================
  // - cuts are defined the same way as bins are defined; be mindful
  //   of what bins you are defining vs. what cuts you are defining.
  //   For example, if you define a Q2 minimum cut, and you also define
  //   Q2 bins below, you may be creating more bins than you actually
  //   need, since the Q2 minimum cut actually defines another bin;
  //   in this case, your Q2 bins effectively define a Q2 minimum.
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


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
