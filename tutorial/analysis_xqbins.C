// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)

/* run in a grid of (x,Q2) 2D bins
 * - various ways to make a grid are demonstrated
 * - observe how the resulting histograms differ in each (x,Q2) bin
 */
void analysis_xqbins(
    TString configFile="tutorial/delphes.config", // input config file
    TString outfilePrefix="tutorial.xqbins"       // output filename prefix
) {

  // setup analysis ========================================
  AnalysisDelphes *A = new AnalysisDelphes(configFile, outfilePrefix);

  //A->maxEvents = 30000; // use this to limit the number of events
  A->SetReconMethod("Ele"); // set reconstruction method
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state


  // define cuts ====================================
  /* - cuts are defined the same way as bins are defined (see below for syntax details and examples);
   *   the `CutDef` class handles bin definitions
   *
   * For example, if you want to apply the y<0.95 cut, do:
   *   A->AddBinScheme("y");
   *   A->BinScheme("y")->BuildBin("Max",0.95);
   * This makes a single y bin, defined by a maximum of 0.95.
   * 
   * If instead you want to make a cut of 0.01<y<0.95, you must do
   *   A->AddBinScheme("y");
   *   A->BinScheme("y")->BuildBin("Range",0.01,0.95);
   * to define a single y bin. If instead you try to do something like
   *   A->AddBinScheme("y");
   *   A->BinScheme("y")->BuildBin("Min",0.01);
   *   A->BinScheme("y")->BuildBin("Max",0.95);
   * you would actually be defining two y bins, one with the minimum and a
   * second with the maximum, which may not be what you want.
   * 
   * You must also be mindful of what other bins you are defining. Suppose you
   * want to apply a Q2>10 GeV2 cut, and then do an analysis in bins of Q2. If
   * you do
   *   A->AddBinScheme("q2");
   *   A->BinScheme("q2")->BuildBin("Min",10); // Q2>10 GeV2 cut
   *   A->BinScheme("q2")->BuildBins( 5, 10,    100,  true ); // 5 bins in range 10-100 GeV, equal width log scale
   * you are actually defining 6 Q2 bins: the 5 you specify with `BuildBins`
   * plus the one you specified with `BuildBin("Min",10)`. In this case, only
   * the third line is needed to apply the Q2>10 GeV2 cut, since your binning
   * range starts at 10 GeV2.
   */ 
  // here are some common cuts we typically use for SIDIS:
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("y");  A->BinScheme("y")->BuildBin("Range",0.01,0.95); // 0.01 < y < 0.95
  A->AddBinScheme("z");  A->BinScheme("z")->BuildBin("Range",0.2,0.9); // 0.2 < z < 0.9
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)



  // set binning scheme ====================================

  // first add what bin schemes you want; you don't have to define
  // bins for all of them, but you need to add at least the ones
  // you plan to use below
  A->AddBinScheme("q2");
  A->AddBinScheme("x");
  A->AddBinScheme("pt");

  // then add the bins that you want
  // - tutorial switch statement: change tutorial number to try out
  //   the different binning implementations
  int tutorialNum = 1;
  switch(tutorialNum) {

    case 0:
      // 1D binning in Q2, equal width in linear scale
      // - arguments of BuildBins are: (numBins, min, max, log-scale bool)
      // - alternatively: BuildBins(TAxis*, log-scale bool)
      // - log-scale is false by default
      A->BinScheme("q2")->BuildBins(  3, 1, 100);
      break;

    case 1:
      // 3x3 grid of (x,Q2) bins, equal width in logarithmic scale
      A->BinScheme("q2")->BuildBins( 3, 1,    100,  true );
      A->BinScheme("x")->BuildBins(  3, 0.01, 1,    true );
      break;

    case 2:
      // alternatively: equal width in linear scale
      A->BinScheme("q2")->BuildBins( 3, 1,    100 );
      A->BinScheme("x")->BuildBins(  3, 0.01, 1   );
      break;

    case 3:
      // custom 2x2 grid (see `CutDef` class for more cut definitions)
      // - arguments of BuildBin are (cutType, a, b), where cutType is
      //   one of the types given in `../src/CutDef.cxx`, and a and b
      //   depend on which cutType
      // - various cutTypes are exemplified here:
      A->BinScheme("q2")->BuildBin("Max",10); // Q2 < 10 GeV2
      A->BinScheme("q2")->BuildBin("Min",10); // Q2 > 10 GeV2
      A->BinScheme("x")->BuildBin("Range", 0.05, 0.2 ); // 0.05 < x < 0.2
      A->BinScheme("x")->BuildBin("CenterDelta", 0.5, 0.1 ); // |x-0.5|<0.1
      break;

    case 4:
      // overlapping Q2 bins, by specifying various Q2 minima
      // - bins are arbitrary and allowed to be overlapping
      A->BinScheme("q2")->BuildBin("Min",1);
      A->BinScheme("q2")->BuildBin("Min",10);
      A->BinScheme("q2")->BuildBin("Min",50);
      break;

    case 5:
      // 3D binning: 2x2 grid of (x,Q2) bins, for two pT bins
      // - you can add more dimensions, but be careful of the curse
      //   of dimensionality
      A->BinScheme("q2")->BuildBins( 2, 1,    100,  true );
      A->BinScheme("x")->BuildBins(  2, 0.01, 1,    true );
      A->BinScheme("pt")->BuildBin("Range", 0.2, 0.5 );
      A->BinScheme("pt")->BuildBin("Range", 0.5, 0.8 );
      break;
  };



  // perform the analysis ==================================
  A->Execute();

  // for reference, here is a print out of HistosDAG
  // - it lists each node, together with its inputs and outputs, which
  //   indicate the connections between the nodes
  //A->GetHistosDAG()->PrintBreadth("HistosDAG Nodes");
};
