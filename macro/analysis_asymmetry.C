R__LOAD_LIBRARY(Largex)

struct WeightsProkudin : public WeightsSivers {
  ProkudinSfSet sf_set;
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

struct WeightsPavia : public WeightsSivers {
  PaviaSfSet sf_set;
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
    TString infiles="datarec/dire_5x41.brian.hiDiv.root", /* delphes tree(s) */
    Double_t eleBeamEn=5, /* electron beam energy [GeV] */
    Double_t ionBeamEn=41, /* ion beam energy [GeV] */
    Double_t crossingAngle=0, /* crossing angle [mrad] */
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
    new WeightsTest()
  });
  A->SetWeights(weights);

  A->writeSimpleTree = true; // write SimpleTree (for one bin)
  //A->maxEvents = 10000; // use this to limit the number of events
  //A->SetReconMethod("Ele"); // recon method (default is "Ele")
  //A->AddFinalState("pipTrack"); // final states (default is "pipTrack" only)

  // perform the analysis ==================================
  A->Execute();
};

