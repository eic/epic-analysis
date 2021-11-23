R__LOAD_LIBRARY(Largex)

struct WeightsTest : public WeightsSivers {
  Double_t Asymmetry(Int_t h, Double_t x, Double_t z, Double_t Q2, Double_t pt) const override {
    return 0.2;
  }
};

void analysis_asymmetry(
    TString infiles="datarec/in.config", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=25, /* crossing angle [mrad] */
    TString outfilePrefix="asymmetry" /* output filename prefix*/
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(
      infiles,
      eleBeamEn,
      ionBeamEn,
      crossingAngle,
      outfilePrefix
      );
  Weights* weights = new WeightsSum({
    new WeightsUniform(),
    new PaviaWeights()
  });
  A->SetWeights(weights);
  A->SetPolT(0.7);
  //A->SetPolL(0.7);
  //A->SetPolB(0.8);

  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  //A->maxEvents = 10000; // use this to limit the number of events
  //A->SetReconMethod("Ele"); // recon method (default is "Ele")
  A->AddFinalState("pipTrack"); // final states (default is "pipTrack" only)

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // perform the analysis ==================================
  A->Execute();
};

