R__LOAD_LIBRARY(Largex)

void analysis_xqbins(
    TString infiles="container/ci/macro/files.config", /* list of input files */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=25, /* crossing angle [mrad] */
    TString outfilePrefix="macro.test" /* output filename prefix*/
) {

  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  // settings
  A->SetReconMethod("Ele");
  A->AddFinalState("pipTrack");

  // cuts
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // bins
  A->AddBinScheme("q2");
  A->AddBinScheme("x");
  A->BinScheme("q2")->BuildBins( 3, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );

  A->Execute();
};
