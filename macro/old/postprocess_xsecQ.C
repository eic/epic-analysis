// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

R__LOAD_LIBRARY(EpicAnalysis)
#include "PostProcessor.h"

// dump counts and average kinematics in Q2 bins
void postprocess_xsecQ(
    TString infile="out/xsecQ.crossCheck.root"
) {

  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  // - specify settings such as diagonalized binnings
  AnalysisDelphes *A = new AnalysisDelphes();
  A->diagonalPtXZ = true;

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(infile);

  // start output file  ========================================
  TString tableFile = P->GetPngDir() + "/table_counts_qbins.txt";
  P->StartTextFile(tableFile,"Counts and averages in bins of Q2");

  // loop over (pt,x,z) bins, diagonalized
  // TODO: when AnalysisDelphes::histSet is generalized, loops like this will
  // be significantly simpler
  for(int bpt : P->GetBinNums("pt")) {
  for(int bx  : P->GetBinNums("x")) {
  for(int bz  : P->GetBinNums("z")) {
    if(A->CheckDiagonal(bpt,bx,bz,-1)) continue; // diagonalize

    // header for this (pt,x,z) bin
    P->AppendToTextFile(tableFile,Form("\nKinematic Bin: %d =========================",bpt));

    // loop over y minima and final states; we will have one table per iteration
    for(int by  : P->GetBinNums("y")) {
    for(int bfs : P->GetBinNums("finalState")) {

      // loop over Q2 bins; these are the rows of the table
      for(int bq : P->GetBinNums("q2")) {

	// skip full bin
	if(P->GetBinCut("q2",bq)->GetCutType()=="Full") continue;

        // ALGORITHM: dump tables of average values, in bins of Q2
        P->DumpAve(
            tableFile,
            A->GetHistosName(bpt,bx,bz,bq,by,bfs),
            "q2");
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
