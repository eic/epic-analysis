// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// make grids of plots of (x,Q2) bins
// - depending on infile, different histograms will be drawn
void postprocess_x_q2(TString infile="out/coverage.fastsim.root") {

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

  // (x,Q2) binning
  auto xBins = P->Op()->GetBinSet("x");
  auto qBins = P->Op()->GetBinSet("q2");
  Int_t numXbins = xBins->GetNumBins();
  Int_t numQbins = qBins->GetNumBins();
  Double_t xMin = xBins->GetMin();
  Double_t xMax = xBins->GetMax();
  Double_t qMin = qBins->GetMin();
  Double_t qMax = qBins->GetMax();

  // 2D array of Histos pointers
  std::vector<std::vector<Histos*>> histosArr(numXbins,std::vector<Histos*>(numQbins));

  // payload operator: find (x,Q2) bin, get (x,Q2) ranges, fill histosArr
  auto fillHistosArr = [&histosArr](NodePath *NP, Histos *H ) {
    Int_t bx = NP->GetBinNode("x")->GetBinNum();
    Int_t bq = NP->GetBinNode("q2")->GetBinNum();
    printf("   bx, bq = %d, %d\n",bx,bq);
    try { histosArr.at(bx).at(bq) = H; }
    catch(const std::out_of_range &e) { cerr << "ERROR: (x,Q2) bin number (" << bx << "," << bq << ") invalid" << endl; };
  };

  // after subloop operator: draw array of plots in (x,Q2) bins
  auto drawHistosArr = [&](NodePath *NP) {
    TString canvName = "xQ2cov_" + NP->BinListName();
    for( TString histName : histList ) {
      P->DrawInBins(
          canvName, histosArr, histName,
          "x", numXbins, xMin, xMax, true,
          "Q^{2}", numQbins, qMin, qMax, true
          );
    };
  };

  // staging and execution
  P->Op()->AfterSubloop( {"x","q2"}, drawHistosArr );
  P->Op()->Payload(fillHistosArr);
  P->Execute();
  P->Finish();
};
