// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)
#include "PostProcessor.h"

// dump counts and average kinematics in pT bins
void postprocess_xsecPT(
    TString infile="out/xsecPT.crossCheck.root"
) {

  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  // - specify settings such as diagonalized binnings
  AnalysisDelphes *A = new AnalysisDelphes();
  A->diagonalXZQ = true;

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // start output file  ========================================
  TString tableFile = P->GetPngDir() + "/table_counts_ptbins.txt";
  P->StartTextFile(tableFile,"Counts and averages in bins of pT");

  // loop over (pt,x,z) bins, diagonalized
  // TODO: when AnalysisDelphes::histSet is generalized, loops like this will
  // be significantly simpler
  for(int bx  : P->GetBinNums("x")) {
  for(int bz  : P->GetBinNums("z")) {
  for(int bq  : P->GetBinNums("q2")) {
    if(A->CheckDiagonal(-1,bx,bz,bq)) continue; // diagonalize

    // header for this (x,z,q2) bin
    P->AppendToTextFile(tableFile,Form("\nKinematic Bin: %d =========================",bx));

    // loop over y minima and final states; we will have one table per iteration
    for(int by  : P->GetBinNums("y")) {
    for(int bfs : P->GetBinNums("finalState")) {

      // loop over pT bins; these are the rows of the table
      for(int bpt : P->GetBinNums("pt")) {

        // skip full bin
        if(P->GetBinCut("pt",bpt)->GetCutType()=="Full") continue;

        // ALGORITHM: dump tables of average values, in bins of pT
        P->DumpAve(
            tableFile,
            A->GetHistosName(bpt,bx,bz,bq,by,bfs),
            "pt");
      };

      // finish ALGORITHM - called after the loop over table rows, so that
      // PostProcessor knows to start a new table for the next y minimum
      P->FinishDumpAve(tableFile);
    }};
  }}};

  // dump final table to stdout
  P->PrintTextFile(tableFile);
  cout << tableFile << " written" << endl;

  // finish
  P->Finish();
};
