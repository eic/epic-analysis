// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Connor Pecar

R__LOAD_LIBRARY(EpicAnalysis)

/*
Script to plot histograms binned in two variables in their respective positions
on those axes. In this script, pT vs Eta plots are plotted in Q2 and 
x bins. This script assumes that the only "Range" bins for Q2 and x are
continuous and of the same size (in regular or log scale), e.g. set 
using the BuildBins function.
*/
  
void postprocess_pTvsEta(
     TString infile = "out/coverage_dagtest.pythia8NCDISS3_10x100_Q21_cross-0.025.root"
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
    P->DrawInBins(canvname, histos_xQ2, "etaVsP", "x", nx, 1e-4, 1, true, "Q^{2}", nq2, .99, 1000, true);
  };

  auto beforefunction = [](){

  };
  P->Op()->Subloop({"x","q2"},beforefunction,drawinxQ2bins);
  P->Op()->Payload(findxQ2bins);
  P->Op()->PrintBreadth("HistosDAG Final Setup");

  P->Execute();
  
  P->Finish();
};
