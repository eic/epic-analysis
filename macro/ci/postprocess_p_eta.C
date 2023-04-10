// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make grids of plots of (p,eta) bins
// - depending on infile, different histograms will be drawn
void postprocess_p_eta(TString infile="out/coverage.fullsim.root") {

  // set histogram lists, based on infile name
  std::vector<TString> histList;
  if(infile.Contains("coverage")) {
    histList.push_back("Q2vsX");
    histList.push_back("etaVsP");
    histList.push_back("phiHvsPhiS");
  }
  else if(infile.Contains("resolution")) {
    histList.push_back("x_Res");
    histList.push_back("y_Res");
    histList.push_back("pT_Res");
    histList.push_back("Q2_Res");
    histList.push_back("Nu_Res");
    histList.push_back("W_Res");
    histList.push_back("phiH_Res");
    histList.push_back("phiS_Res");
    histList.push_back("z_Res");
    histList.push_back("mX_Res");
    histList.push_back("xF_Res");
  } else {
    fprintf(stderr,"ERROR: histList not defined for specified infile name\n");
    return;
  };

  // build DAG
  PostProcessor *P = new PostProcessor(infile);

  // (p,eta) binning
  auto pBins     = P->Op()->GetBinSet("p");
  auto eBins     = P->Op()->GetBinSet("eta");
  Int_t numPbins = pBins->GetNumBins();
  Int_t numEbins = eBins->GetNumBins();
  Double_t pMin  = pBins->GetMin();
  Double_t pMax  = pBins->GetMax();
  Double_t eMin  = eBins->GetMin();
  Double_t eMax  = eBins->GetMax();

  // 2D array of Histos pointers
  std::vector<std::vector<Histos*>> histosArr(numPbins,std::vector<Histos*>(numEbins));

  // payload operator: find (p,eta) bin, get (p,eta) ranges, fill histosArr
  auto fillHistosArr = [&histosArr](NodePath *NP, Histos *H ) {
    Int_t bp = NP->GetBinNode("p")->GetBinNum();
    Int_t be = NP->GetBinNode("eta")->GetBinNum();
    printf("   bp, be = %d, %d\n",bp,be);
    try { histosArr.at(bp).at(be) = H; }
    catch(const std::out_of_range &e) { cerr << "ERROR: (p,eta) bin number (" << bp << "," << be << ") invalid" << endl; };
  };

  // after subloop operator: draw array of plots in (p,eta) bins
  auto drawHistosArr = [&](NodePath *NP) {
    TString canvName = "pEtaCov_" + NP->BinListName();
    for( TString histName : histList ) {
      P->DrawInBins(
          canvName, histosArr, histName,
          "p",      numPbins, pMin, pMax, true,
          "#eta",   numEbins, eMin, eMax, false
          );
    };
  };

  // staging and execution
  P->Op()->AfterSubloop( {"p","eta"}, drawHistosArr );
  P->Op()->Payload(fillHistosArr);
  P->Execute();
  P->Finish();
};
