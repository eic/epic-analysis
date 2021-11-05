R__LOAD_LIBRARY(Largex)

struct WeightsSfSet : public WeightsSivers {
  SfSet const& sf_set;
  WeightsSfSet(SfSet const& sf_set) : sf_set(sf_set) { }
  Double_t Asymmetry(Double_t x, Double_t z, Double_t Q2, Double_t pt) const override {
    Double_t fuu = sf_set.F_UUT(Hadron::PI_P, x, z, Q2, pt*pt);
    Double_t fut = sf_set.F_UTT_sin_phih_m_phis(Hadron::PI_P, x, z, Q2, pt*pt);
    if (!TMath::Finite(fuu) || !TMath::Finite(fut)) {
        fuu = 1.;
        fut = 0.;
    }
    if (TMath::Abs(fut / fuu) > 1.) {
        fuu = 1.;
        fut = TMath::Sign(1., fut / fuu);
    }
    return fut / fuu;
  }
};

struct WeightsTest : public WeightsSivers {
  Double_t Asymmetry(Double_t x, Double_t z, Double_t Q2, Double_t pt) const override {
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

  //ProkudinSfSet sf_set;
  //PaviaSfSet sf_set;
  TorinoSfSet sf_set;
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
    new WeightsSfSet(sf_set)
  });
  A->SetWeights(weights);

  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  //A->maxEvents = 10000; // use this to limit the number of events
  //A->SetReconMethod("Ele"); // recon method (default is "Ele")
  //A->AddFinalState("pipTrack"); // final states (default is "pipTrack" only)

  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // perform the analysis ==================================
  A->Execute();
};

