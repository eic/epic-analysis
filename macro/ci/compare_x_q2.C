R__LOAD_LIBRARY(Largex)

// make grids of plots of (x,Q2) bins, comparing data from the infiles
// - depending on infile, different histograms will be drawn
void compare_x_q2(
    TString infile0="out/resolution.fullsim.Ele.root",
    TString infile1="out/resolution.fullsim.DA.root"
    ) {

  // histograms ================================================================
  // - set histogram lists, based on infile name
  std::vector<TString> histList;
  if(infile0.Contains("coverage")) {
    histList.push_back("Q2vsX");
    histList.push_back("x");
    histList.push_back("y");
    histList.push_back("W");
    histList.push_back("pLab");
    histList.push_back("pTlab");
    histList.push_back("etaLab");
    histList.push_back("phiLab");
    histList.push_back("z");
    histList.push_back("pT");
    histList.push_back("qT");
    histList.push_back("qTq");
    histList.push_back("mX");
    histList.push_back("phiH");
    histList.push_back("phiS");
    histList.push_back("phiHvsPhiS");
    histList.push_back("etaVsP");
  }
  else if(infile0.Contains("resolution")) {
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


  // setup ================================================================

  std::vector<TString> infiles;
  infiles.push_back(infile0);
  infiles.push_back(infile1);
  bool first=true;
  Int_t numXbins, numQbins;
  Double_t xMin, xMax, qMin, qMax;

  PostProcessor *P0; // first dag
  std::vector<PostProcessor*> procs; // additional dags
  PostProcessor *P; // pointer

  // get dags and binning
  for(auto infile : infiles) {
    P = new PostProcessor(infile);
    auto xBins = P->Op()->GetBinSet("x");
    auto qBins = P->Op()->GetBinSet("q2");
    if(first) {
      first = false;
      P0 = P;
      numXbins = xBins->GetNumBins();
      numQbins = qBins->GetNumBins();
      xMin = xBins->GetMin();
      xMax = xBins->GetMax();
      qMin = qBins->GetMin();
      qMax = qBins->GetMax();
    } else {
      procs.push_back(P);
      if(numXbins != xBins->GetNumBins() || numQbins != qBins->GetNumBins()) {
        cerr << "ERROR: files have differing bins" << endl;
        return;
      }
    }
  }

  // list of 2D arrays of Histos pointers; each element of the list will be compared
  int numFiles = infiles.size();
  std::vector<std::vector<std::vector<Histos*>>> histosArrList(
      numFiles,
      std::vector<std::vector<Histos*>>(
        numXbins,
        std::vector<Histos*>(numQbins)
        )
      );


  // operators ================================================================
  // - for P0 dag

  // payload: find (x,Q2) bin, and insert into histosArrList, for each infile
  auto fillHistosArr = [&histosArrList,procs](NodePath *NP, Histos *H ) {
    Int_t bx = NP->GetBinNode("x")->GetBinNum();
    Int_t bq = NP->GetBinNode("q2")->GetBinNum();
    Int_t pc=0;
    printf("   bx, bq = %d, %d\n",bx,bq);
    try { 
      histosArrList.at(0).at(bx).at(bq) = H; // first, insert the Histos* of P0
      for(auto proc : procs) { // then insert the Histos* of each dag in procs
        histosArrList.at(++pc).at(bx).at(bq) = proc->Op()->GetHistosExternal(NP);
      }
    }
    catch(const std::out_of_range &e) { 
      fprintf(stderr,"ERROR: invalid bin number (file,x,Q2) = (%d,%d,%d)\n",pc,bx,bq);
    }
  };

  // after subloop operator: draw array of plots in (x,Q2) bins
  auto drawHistosArr = [&](NodePath *NP) {
    TString canvName = "xQ2cov_" + NP->BinListName();
    for( TString histName : histList ) {
      P0->DrawInBins(
          canvName, histosArrList, histName,
          "x",      numXbins,      xMin,     xMax, true,
          "Q^{2}",  numQbins,      qMin,     qMax, true
          );
    };
  };

  // staging and execution
  P0->Op()->AfterSubloop( {"x","q2"}, drawHistosArr );
  P0->Op()->Payload(fillHistosArr);
  P0->Execute();
  P0->Finish();
};
