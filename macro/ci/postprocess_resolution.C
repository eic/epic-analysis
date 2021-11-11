R__LOAD_LIBRARY(Largex)

// make resolution plots
void postprocess_resolution(
    TString infile="out/resolution.root"
){

  // cleanup old image files
  gROOT->ProcessLine(".! rm -v out/resolution.images/*.png");
  gROOT->ProcessLine(".! rm -v out/resolution.images/*.pdf");

  // build DAG
  PostProcessor *P = new PostProcessor(infile);
  //P->Op()->PrintBreadth("HistosDAG Initial Setup");

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
  auto drawHistosArr = [&histosArr, &P, &numXbins, &numQbins, &xMin, &xMax, &qMin, &qMax](NodePath *NP) {
    TString canvName = "xQ2cov_" + NP->BinListName();
    for( TString histName : {"x_Res","y_Res","pT_Res","Q2_Res","phiH_Res","phiS_Res","phiHvsPhiS"} ) {
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
