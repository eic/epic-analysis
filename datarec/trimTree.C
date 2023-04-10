// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// copy numEntries entries of a TTree to another file
void trimTree(TString inFileName, TString treeName, Long64_t numEntries) {
  TFile *inFile = new TFile(inFileName,"READ");
  TString outFileName = inFileName;
  outFileName(TRegexp("root$")) = "trimmed.root";
  cout << "outFileName = " << outFileName << endl;
  TFile *outFile = new TFile(outFileName,"RECREATE");
  TTree *inTree = (TTree*) inFile->Get(treeName);
  TTree *outTree = inTree->CloneTree(numEntries);
  outTree->Write();
  outFile->Close();
  inFile->Close();
};
