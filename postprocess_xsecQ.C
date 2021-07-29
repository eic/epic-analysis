R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// dump counts and average kinematics in Q bins
void postprocess_xsecQ() {

  // instantiate empty analysis object ================================
  // - needed for some general information about binning
  // - specify settings such as diagonalized binnings
  Analysis *A = new Analysis();
  A->diagonalPtXZ = true;

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(
      "out/histos.dire_5x41.brian.hiAcc.10milEvents.root", /* histograms file */
      false
      );

  // start output file  ========================================
  TString tableFile = P->GetPngDir() + "/table_counts_qbins.txt";
  P->StartTextFile(tableFile,"Counts and averages in bins of Q");

  // loop over (pt,x,z) bins, diagonalized
  for(int bpt : P->GetBinNums("pt")) {
  for(int bx  : P->GetBinNums("x")) {
  for(int bz  : P->GetBinNums("z")) {
    if(A->CheckDiagonal(bpt,bx,bz,-1)) continue; // diagonalize

    // header for this (pt,x,z) bin
    P->AppendToTextFile(tableFile,Form("\nKinematic Bin: %d =========================",bpt));

    // loop over y minima and final states; we will have one table per iteration
    for(int by  : P->GetBinNums("y")) {
    for(int bfs : P->GetBinNums("finalState")) {

      // loop over Q bins; these are the rows of the table
      for(int bq : P->GetBinNums("q")) {

        // ALGORITHM: dump tables of average values, in bins of "Q"
        P->DumpAve(
            tableFile,
            A->GetHistosName(bpt,bx,bz,bq,by,bfs),
            "q");
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


  // TODO: make these additional postprocessor macros

  /*
  // ratios y>y_min / y>0
  TString kinBinStr = "histos_pipTrack_pt0_x0_z0_q0";
  ResetVars();
  // draw
  for(int y=1; y<=3; y++) {
    DrawRatios(
        Form("yRatio%d",y),
        kinBinStr+Form("_y%d",y),
        kinBinStr+"_y0"
        );
  };
  // write summary
  outfile->cd("/");
  outfile->mkdir("summary_yRatio");
  outfile->cd("summary_yRatio");
  for(auto const& kv : summaryCanvMap) {
    kv.second->Write();
    kv.second->Print(pngDir+"/summary_yRatio_"+kv.first+".png");
  };
  outfile->cd("/");
  */

  // -------------------------------------------------

  /*
  // dump counts in Q bins
  ResetVars();
  Int_t ycut=0; // choose y cut
  TString setStr,histStr;
  TString yieldOutput = pngDir + "/table_counts.txt";
  gSystem->RedirectOutput(yieldOutput,"w");
  cout << "Counts in bins of Q" << endl;
  gSystem->RedirectOutput(0);
  for(int b=0; b<=3; b++) {
    setStr = Form("histos_pipTrack_pt%d_x%d_z%d_q%d_y%d",b,b,b,0,ycut);
    DrawSingle(Form("xsec_%d",b),setStr,"Q_xsec");
    DumpHist(yieldOutput,setStr,"Q");
  };
  // write summary
  outfile->cd("/");
  summaryCanv->Write();
  summaryCanv->Print(pngDir+"/summary_xsec.png");
  // dump
  gROOT->ProcessLine(".! cat "+yieldOutput);
  cout << yieldOutput << " written" << endl;
  */
      
  // -------------------------------------------------


};
