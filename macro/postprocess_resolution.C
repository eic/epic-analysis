R__LOAD_LIBRARY(Largex)

// make resolution plots
// - adapted from `postprocess_pTvsEta.C`
void postprocess_resolution(
    TString infile="out/resolution.root"
){
  
  PostProcessor *P = new PostProcessor(infile);
  P->Op()->PrintBreadth("HistosDAG Initial Setup");

  // number of bins in x and Q2
  int nx = P->Op()->GetBinSet("x")->GetNumBins();
  int nq2 = P->Op()->GetBinSet("q2")->GetNumBins();

  // just counters for filling Histos vector
  int xbin = 0;
  int q2bin = 0;
  
  // initialize this 2D vector to be some large size
  std::vector<std::vector<Histos*>> histos_xQ2(30,std::vector<Histos*>(30));
  
  auto findxQ2bins = [&histos_xQ2,&P,&xbin,&q2bin,nx,nq2](Histos *H ){
    histos_xQ2[xbin][q2bin] = H;    
    q2bin++;
    if(q2bin == nq2){
      q2bin=0; xbin++; 
      if(xbin == nx) xbin = 0;
    }
  };
  
  auto drawinxQ2bins = [&histos_xQ2, &P, &nx, &nq2](NodePath *bins){    
    TString canvname = "xQ2cov_"; //+ bins->GetVar
    for(Node *bin: bins->GetBinNodes()){
      if(bin->GetVarName() == "finalState"){
        canvname+=bin->GetID();
        canvname+="_";
      }
      if(bin->GetVarName() == "z"){
        canvname+=bin->GetID();
        canvname+="_";
      }
    }

    double xMin = 1e-4;
    double xMax = 1;
    double q2Min = 0.99;
    double q2Max = 1000;

    // loop over resolution histograms (see ../src/Analysis.cxx `DefineHist*` calls 
    // for available histograms, or add your own there)
    for( TString histname : {"x_Res","y_Res","Q2_Res","phiH_Res","phiS_Res","phiHvsPhiS"} ) {
      P->DrawInBins(
          canvname, histos_xQ2, histname,
          "x", nx, xMin, xMax, true,
          "Q^{2}", nq2, q2Min, q2Max, true
          );
    };
  };

  auto beforefunction = [](){
  };

  P->Op()->Subloop({"x","q2"},beforefunction,drawinxQ2bins);
  P->Op()->Payload(findxQ2bins);
  //P->Op()->PrintBreadth("HistosDAG Final Setup");

  P->Execute();
  
  P->Finish();
};
