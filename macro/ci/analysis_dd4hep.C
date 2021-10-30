R__LOAD_LIBRARY(Largex)

// full simulation (dd4hep) test
void analysis_dd4hep(
    TString infiles="macro/ci/s3files.list", /* list of input files (S3 URLs, plus other columns) */ // TODO: UPDATE s3 LIST
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=25, /* crossing angle [mrad] */
    TString outfilePrefix="fullsim.coverage" /* output filename prefix*/
) {

  // setup analysis ========================================
  // - define `AnalysisDD4hep` instead of `AnalysisDelphes`
  AnalysisDD4hep *A = new AnalysisDD4hep(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  //A->maxEvents = 300000; // use this to limit the number of events
  A->writeSimpleTree = true;

  // Set scatt. electron cuts
  A->SetEleEnergyThreshold(eleBeamEn * 0.1);  // default is 10% of beamE
  A->SetIsoCut(0.1);        // default is 10%
  A->SetIsoConeRadius(1.0); // default is 1.0

  // set reconstruction method and final states =============================
  // - see `Analysis` constructor for methods (or other tutorials)
  A->SetReconMethod("Ele");
  A->AddFinalState("pipTrack");

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  A->AddBinScheme("q2");
  A->AddBinScheme("x");

  // 3x2 grid of (x,Q2) bins, equal width in logarithmic scale
  A->BinScheme("q2")->BuildBins( 3, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );

  // perform the analysis ==================================
  A->Execute();
}
