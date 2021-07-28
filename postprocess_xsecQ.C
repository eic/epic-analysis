R__LOAD_LIBRARY(Largex)
#include "PostProcessor.h"

// dump counts and average kinematics in Q bins
void postprocess_xsecQ() {

  // setup postprocessor ========================================
  PostProcessor *P = new PostProcessor(
      "out/histos.dire_5x41.brian.hiDiv.root", /* histograms file */
      false
      );


  // start output file  ========================================
  TString setStr,histStr;
  TString aveOutput = P->GetPngDir() + "/table_counts_qbins.txt";
  gSystem->RedirectOutput(aveOutput,"w");
  cout << "Counts and averages in bins of Q" << endl;
  gSystem->RedirectOutput(0);

  // loop over relevant bins
  // TODO: `BinSet` objects are now streamed to root files; use their iterators instead
  //
  for(int b=0; b<=2; b++) {
    gSystem->RedirectOutput(aveOutput,"a");
    cout << endl << "Kinematic Bin: " << b 
         << " ============================================== " << endl;
    gSystem->RedirectOutput(0);
    for(int y=0; y<=3; y++) {

      P->ResetVars(); // starting a new table, call this to reprint headers

      // loop over Q bins
      for(int q=1; q<=10; q++) {
        setStr = Form("histos_pipTrack_pt%d_x%d_z%d_q%d_y%d",b,b,b,q,y);
        P->DumpAve(aveOutput,setStr,"q");
      };

      // format the table (stored in `tmp` file) and stream to main output file
      P->Columnify(aveOutput+".tmp",aveOutput);

      // remove `tmp` file
      gROOT->ProcessLine(".! rm "+aveOutput+".tmp");
    };
  };

  // print dump to stdout
  gROOT->ProcessLine(".! cat "+aveOutput);
  cout << aveOutput << " written" << endl;

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
