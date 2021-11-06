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

  // get number of bins in x and Q2
  Int_t numXbins = P->Op()->GetBinSet("x")->GetNumBins();
  Int_t numQbins = P->Op()->GetBinSet("q2")->GetNumBins();
  Double_t xMin = 1;
  Double_t xMax = 0;
  Double_t qMin = 1e6;
  Double_t qMax = 0;
  
  // 2D array of Histos pointers
  std::vector<std::vector<Histos*>> histosArr(numXbins,std::vector<Histos*>(numQbins));
  
  // payload operator: find (x,Q2) bin, get (x,Q2) ranges, fill histosArr
  auto fillHistosArr = [&histosArr,&xMin,&xMax,&qMin,&qMax](NodePath *NP, Histos *H ) {
    auto xBin = NP->GetBinNode("x");
    auto qBin = NP->GetBinNode("q2");
    xMin = xBin->GetCut()->GetMin() < xMin ? xBin->GetCut()->GetMin() : xMin;
    xMax = xBin->GetCut()->GetMax() > xMax ? xBin->GetCut()->GetMax() : xMax;
    qMin = qBin->GetCut()->GetMin() < qMin ? qBin->GetCut()->GetMin() : qMin;
    qMax = qBin->GetCut()->GetMax() > qMax ? qBin->GetCut()->GetMax() : qMax;
    Int_t bx = xBin->GetBinNum();
    Int_t bq = qBin->GetBinNum();
    try { histosArr.at(bx).at(bq) = H; }
    catch(const std::out_of_range &e) { cerr << "ERROR: (x,Q2) bin number (" << bx << "," << bq << ") invalid" << endl; };
  };

  // after subloop operator: draw array of plots in (x,Q2) bins
  auto drawHistosArr = [&histosArr, &P, &numXbins, &numQbins](NodePath *NP) { // after subloop

    TString canvName = "xQ2cov_" + NP->BinListName();
    cout << "canvName = " << canvName << endl;

    // loop over resolution histograms (see ../src/Analysis.cxx `DefineHist*` calls 
    // for available histograms, or add your own there)
    for( TString histName : {"x_Res","y_Res","pT_Res","Q2_Res","phiH_Res","phiS_Res","phiHvsPhiS"} ) {
      P->DrawInBins(
          canvName, histosArr, histName,
          "x", numXbins, xMin, xMax, true,
          "Q^{2}", numQbins, qMin, qMax, true
          );
    };
  };

  P->Op()->Payload(fillHistosArr); 
  P->Op()->AfterSubloop( {"x","q2"}, drawHistosArr );

  P->Execute();
  P->Finish();
};
