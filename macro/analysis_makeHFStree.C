// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Connor Pecar

R__LOAD_LIBRARY(EpicAnalysis)

// Produce HFSTree: tree containing all information for kinematic reconstruction
// studies. Includes:
//    -true beam and scattered electron information
//    -reco. scattered electron
//    -std::vector of momentum, energy, and PID for each hadronic final state (HFS) particle
//    -std::vector of momentum, energy, and PID for select final state
//     tracks for SIDIS kinematics (e.g. here pipTrack)
// Does not place any cuts besides requiring a reconstructed scattered electron
// and at least 1 reconstructed HFS particle
void analysis_makeHFStree(
    TString configFile="datarec/in.config", /* delphes tree(s) */
    TString outfilePrefix="resolutions" /* output filename prefix*/
) {

  //outfilePrefix+="_DA";
  // setup analysis ========================================
  AnalysisEpic *A = new AnalysisEpic(
      configFile,
      outfilePrefix
      );
  
  //A->maxEvents = 1000000; // use this to limit the number of events
  A->SetReconMethod("ele"); // set reconstruction method

  A->writeHFSTree = true; // write HFSTree (for one bin)
  A->AddFinalState("pipTrack"); // pion final state
  //A->AddFinalState("KpTrack"); // kaon final state


  // define cuts ====================================
  A->AddBinScheme("w");  A->BinScheme("w")->BuildBin("Min",3.0); // W > 3 GeV
  A->AddBinScheme("xF"); A->BinScheme("xF")->BuildBin("Min",0.0); // xF > 0
  A->AddBinScheme("ptLab");  A->BinScheme("ptLab")->BuildBin("Min",0.1); // pT_lab > 0.1 GeV (tracking limit)


  // set binning scheme ====================================
  // z ranges
  A->AddBinScheme("z");
  A->BinScheme("z")->BuildBin("Min",0.2); // needed?

  // y minima
  A->AddBinScheme("y");
  A->BinScheme("y")->BuildBin("Full"); // a bin with no y-cut

  // perform the analysis ==================================
  A->Execute();

};
