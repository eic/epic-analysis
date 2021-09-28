R__LOAD_LIBRARY(Largex)

void analysis_xqbins(
    TString infiles="/data/ciTest_5x41.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
    TString outfilePrefix="macro.test" /* output filename prefix*/
) {

  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );

  A->SetReconMethod("Ele");
  A->AddFinalState("pipTrack");

  A->AddBinScheme("q2");
  A->AddBinScheme("x");
  A->BinScheme("q2")->BuildBins( 3, 1,    100,  true );
  A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );

  A->Execute();
};
