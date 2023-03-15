// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

// analysis in bins of (x,Q2); tests modifying bin boundaries of Histos histograms
// according to the binning scheme
void analysis_modify_histos(
    TString configFile,
    TString outfilePrefix,
    TString dataSource="epic",
    TString reconMethod="Ele"
    )
{

  // setup analysis ========================================
  Analysis *A;
  if   (dataSource.Contains("epic"))     A = new AnalysisEpic(   configFile, outfilePrefix );
  else if(dataSource.Contains("athena")) A = new AnalysisAthena( configFile, outfilePrefix );
  else if(dataSource.Contains("ecce"))   A = new AnalysisEcce(   configFile, outfilePrefix );
#ifndef EXCLUDE_DELPHES
  else A = new AnalysisDelphes( configFile, outfilePrefix );
#endif

  A->SetReconMethod(reconMethod); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  // A->writeSimpleTree = true;

  // define cuts ===========================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)

  // set binning scheme ====================================
  A->AddBinScheme("x");  A->BinScheme("x")->BuildBins(  6, 0.001, 1,    true );
  A->AddBinScheme("q2"); A->BinScheme("q2")->BuildBins( 4, 1,     3000, true );


  //////////////////////////////////////////////////////////
  // modify histograms according to the binning scheme

  // get the HistosDAG (adage)
  A->Prepare(); // (to initialize the DAG with the above bins)
  auto HD = A->GetHistosDAG();

  // define a function to modify histogram named `hist_name` to
  // have the bin boundaries of the bin named `var_name`; this will
  // run on every possible multi-dimensional bin
  auto mod_bounds = [&HD] (TString hist_name, TString var_name) {

    // payload function, to be run on every multi-dim. bin
    auto payload = [&hist_name,&var_name] (NodePath *NP, Histos *H ) {

      // get the histogram and bin
      auto hist_orig   = H->Hist(hist_name);
      auto hist_config = H->GetHistConfig(hist_name);
      auto var_bin     = NP->GetBinNode(var_name);
      if(hist_orig==nullptr or hist_config==nullptr or var_bin==nullptr)
        return;

      // get the binning info // FIXME: only works for 1D
      auto n_bins  = hist_orig->GetXaxis()->GetNbins();
      auto var_min = var_bin->GetCut()->GetMin();
      auto var_max = var_bin->GetCut()->GetMax();

      // get histogram info
      auto hist_orig_name  = hist_orig->GetName();
      auto hist_orig_title = hist_orig->GetTitle();

      // create new histogram, replacing the old one
      decltype(hist_orig) hist_new;
      switch(hist_orig->GetDimension()) {
        case 1:
          // FIXME: `Warning in <TFile::Append>: Replacing existing TH1`
          hist_new = new TH1D(hist_orig_name, hist_orig_title, n_bins, var_min, var_max);
          if(hist_config->logx) BinSet::BinLog(hist_new->GetXaxis()); // preserve equal bin-width in log scale setting
          H->ReplaceHist(hist_name, hist_new);
          break;
        default:
          fmt::print(stderr,"ERROR: {}-dimensional histogram not fully supported yet\n");
      }
    };

    // run the payload on every multi-dim. bin
    HD->Payload(payload);
    HD->ExecuteAndClearOps();
  };

  // make 1D x distribution have a range that matches the x-bin boundaries
  mod_bounds("x", "x");

  // repeat for Q2; note the variable and histogram names differ
  mod_bounds("Q2", "q2");

  //////////////////////////////////////////////////////////


  // perform the analysis ==================================
  A->Execute();
};
