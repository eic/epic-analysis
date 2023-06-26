// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Connor Pecar, Duane Byer

#include "PostProcessor.h"

ClassImp(PostProcessor)

using std::cout;
using std::cerr;
using std::endl;

// constructor
PostProcessor::PostProcessor(
  TString infileN_,
  TString outfileN_
)
  : infileN(infileN_)
  , outfileN(outfileN_)
{

  // settings
  // - summary canvas formatting
  summaryColor[0] = kRed;     summaryStyle[0] = kFullCircle;
  summaryColor[1] = kGreen+1; summaryStyle[1] = kFullTriangleUp;
  summaryColor[2] = kBlue;    summaryStyle[2] = kFullTriangleDown;
  // - general
  gStyle->SetOptStat(0);
  cout << std::fixed; // fixed float precision
  cout << std::setprecision(3);

  // set up input and output files
  cout << "-- postprocess histograms from " << infileN << endl;
  if(outfileN=="") outfileN = infileN; // default
  outfileN(TRegexp("\\.root$")) = "";
  pngDir = outfileN+".images";
  outfileN += ".canvas.root";
  outfile = new TFile(outfileN,"RECREATE");
  infile = new TFile(infileN,"READ");
  gROOT->ProcessLine(".! mkdir -p "+pngDir);

  // rebuild DAGs
  HD = new HistosDAG();
  HD->BuildFromFile(infile);

  // initialize algorithm-specific vars
  this->ResetVars();

};


//=========================================================================
// cleanup and close open files and streams
void PostProcessor::Finish() {
  infile->Close();
  outfile->Close();
  cout << outfileN << " written." << endl;
  cout << pngDir << "/ images written." << endl;
};


//=========================================================================
// reset global variables, such as summary canvases
void PostProcessor::ResetVars() {
  nsum=0;
  ndump=0;
  summaryCanvMap.clear();
  varList.clear();
};


//=========================================================================
/* ALGORITHM: dump histogram counts
 * - output will be appended to `countsFiles`
 * - so far only implemented for 1D
 */
void PostProcessor::DumpHist(TString datFile, TString histSet, TString varName) {
  cout << "dump " << histSet << " : " << varName << " to " << datFile << endl;
  Histos *H = (Histos*) infile->Get(histSet);
  TH1 *hist = H->Hist(varName);
  if(hist->GetDimension()>1) return;
  TString histTformatted = hist->GetTitle();
  histTformatted.ReplaceAll("#in"," in ");
  histTformatted.ReplaceAll("#pm","+-");

  gSystem->RedirectOutput(datFile,"a");
  cout << endl;
  cout << "Histogram: " << histTformatted << endl;

  gSystem->RedirectOutput("tempo","w");
  cout << "bin"
       << " " << varName << "_min"
       << " " << varName << "_max"
       << " " << "counts"
       << " " << "error"
       << " " << endl;
  for(int b=1; b<=hist->GetNbinsX(); b++) {
    cout << b
         << " " << hist->GetBinLowEdge(b)
         << " " << hist->GetBinLowEdge(b+1)
         << " " << hist->GetBinContent(b)
         << " " << hist->GetBinError(b)
         << endl;
  };

  gSystem->RedirectOutput(datFile,"a");
  gROOT->ProcessLine(".! cat tempo | column -t");
  gSystem->RedirectOutput(0);
  gROOT->ProcessLine(".! rm tempo");

};

//=========================================================================

/* ALGORITHM: dump averages of all histograms from `histSet` to the file `datFile`
 * - it is best to call this function in a loop, this function will only dump
 *   one line of information; if it's the first time you called it, it will
 *   dump header information as well
 * - use `cutName` to specify a cut defintion you want to print; e.g., loop
 *   over x bins and set the `cutName` to the x cut, to include columns for x
 *   bin boundaries
 * - when done looping, call `FinishDumpAve(datFile)`
 */
void PostProcessor::DumpAve(TString datFile, Histos *H, TString cutName) {
  cout << "dump averages from " << H->GetSetName()
       << " to " << datFile << endl;

  // get associated cut, and print header info if it's the first time
  dumpCut = nullptr;
  if(ndump==0) {
    gSystem->RedirectOutput(datFile,"a");
    cout << endl << "Histogram Set:";
  };
  for(CutDef *cut : H->CutDefList) {
    if(cutName.CompareTo(cut->GetVarName(),TString::kIgnoreCase)==0) {
      dumpCut = cut;
    }
    else if(ndump==0) {
      TString cutT = cut->GetCutTitle();
      cutT.ReplaceAll("#in"," in ");
      cutT.ReplaceAll("#pm","+-");
      cout << "   " + cutT;
    };
  };
  if(ndump==0) {
    cout << endl;
    cout << "Averages in bins of " << dumpCut->GetVarTitle() << ":" << endl;
    cout << "------------------------" << endl;
    gSystem->RedirectOutput(0);
  };
  if(dumpCut==nullptr) {
    cerr << "ERROR: in DumpAve, cannot find cut def" << endl;
    return;
  };

  // start temporary file for output; this is so you can use `Columnify` later
  // for pretty print to full `datFile`
  gSystem->RedirectOutput(datFile+".tmp",ndump>0?"a":"w");

  // if this is the first time, populate ordered list of histograms; we keep a
  // local list, in case you call this function on a `Histos` with different
  // ordering; we want to maintain order for printout
  if(ndump==0) {
    // header for bin columns
    cout << "bin"
         << " " << dumpCut->GetVarTitle() << "_min"
         << " " << dumpCut->GetVarTitle() << "_max"
         << " counts";
    // loop over 1D histograms, build `varList`,
    // and print header for ave columns
    for(TString varName : H->VarNameList) {
      if(H->Hist(varName)->GetDimension()==1) {
        if(varName.Contains("xsec",TString::kIgnoreCase)) continue;
        if(varName.Contains("Res",TString::kIgnoreCase)) continue;
        if(varName.Contains("RvG",TString::kIgnoreCase)) continue;
        if(H->Hist(varName)->GetEntries()==0) continue;
        varList.push_back(varName);
        cout << " <" + varName + ">";
      };
    };
    // header newline
    cout << endl;
  };

  // loop over local list of variables, and print their averages
  TH1 *hist;
  cout << ndump
       << " " << dumpCut->GetMin()
       << " " << dumpCut->GetMax();
  Bool_t first=true;
  Double_t counts;
  for(TString varName : varList) {
    hist = H->Hist(varName);
    if(first) {
      counts = hist->GetEntries();
      cout << " " << counts;
      first=false;
    }
    else if(hist->GetEntries()!=counts) {
      cerr << "WARNING: mismatch counts in " << varName << " histogram" << endl;
    };
    cout << " " << hist->GetMean();
  };
  cout << endl;
  gSystem->RedirectOutput(0);


  // iterate row counter
  ndump++;

};
void PostProcessor::FinishDumpAve(TString datFile) {
  this->Columnify(datFile+".tmp",datFile);
  gROOT->ProcessLine(".! rm "+datFile+".tmp");
  this->ResetVars();
};


//=========================================================================
/* ALGORITHM: draw a single histogram to a canvas, and write it
 * - `histName` is the name of the histogram given in Histos
 * - `drawFormat` is the formatting string passed to TH1::Draw()
 * - `profileAxis` will draw a TProfile on the specified axis, for 2D dists
 *   - 0=disabled(default), 1=x-axis, 2=y-axis
 * - `profileOnly` draw only the TProfile, for 2D dists
 */
void PostProcessor::DrawSingle(Histos *H, TString histName, TString drawFormat, Int_t profileAxis, Bool_t profileOnly) {
  cout << "draw single plot " << histName << "..." << endl;
  TH1 *hist = H->Hist(histName);
  if(hist==nullptr) {
    cerr << "ERROR: cannot find histogram " << histName << endl;
    return;
  };
  TString canvN = "canv_"+histName+"___"+H->GetSetName();
  TCanvas *canv = new TCanvas(canvN,canvN,dimx,dimy);
  hist->Draw(drawFormat);
  if(hist->GetMinimum()>=0 && hist->GetDimension()==1 && H->GetHistConfig(histName)->logy==false)
    hist->GetYaxis()->SetRangeUser(0,hist->GetMaximum()*1.1); // do not suppress zero

  canv->SetGrid(1,1);
  canv->SetLogx(H->GetHistConfig(histName)->logx);
  canv->SetLogy(H->GetHistConfig(histName)->logy);
  canv->SetLogz(H->GetHistConfig(histName)->logz);
  canv->SetBottomMargin(0.15);
  canv->SetLeftMargin(0.15);
  canv->SetRightMargin(0.15);

  // profile for 2D plot
  if(profileAxis>0 && hist->GetDimension()==2) {
    TProfile *prof;
    switch(profileAxis) {
      case 1: prof = ((TH2*)hist)->ProfileX(); break;
      case 2: prof = ((TH2*)hist)->ProfileY(); break;
    };
    prof->SetLineColor(kBlack);
    prof->SetLineWidth(3);
    if(profileOnly) {
      prof->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
      prof->GetYaxis()->SetRangeUser(hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
      prof->Draw();
    } else prof->Draw("same");
    outfile->cd("/");
    prof->Write();
  }

  // profile for 3D plot
  if(profileAxis>0 && hist->GetDimension()==3) {
    TProfile2D *prof = ((TH3*)hist)->Project3DProfile("yx"); // FIXME: generalize for other axes
    prof->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
    prof->GetYaxis()->SetRangeUser(hist->GetYaxis()->GetXmin(),hist->GetYaxis()->GetXmax());
    prof->GetZaxis()->SetRangeUser(hist->GetZaxis()->GetXmin(),hist->GetZaxis()->GetXmax());
    prof->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
    prof->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
    prof->GetZaxis()->SetTitle(hist->GetZaxis()->GetTitle());
    prof->Draw(drawFormat);
    outfile->cd("/");
    prof->Write();
  }

  canv->Print(pngDir+"/"+canvN+".png");
  outfile->cd("/");
  canv->Write();
};

// OLD VERSION: 
void PostProcessor::DrawSingle(TString histSet, TString histName) {

  Histos *H = (Histos*) infile->Get(histSet);
  TH1 *hist = H->Hist(histName,true);
  Hist4D *hist4 = H->Hist4(histName,true);

  if (hist != nullptr) {
    TString canvN = "canv_"+histName+"_"+H->GetSetName();
    TCanvas *canv = new TCanvas(canvN,canvN, dimx, dimy);

    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(1.0);
    hist->SetLineWidth(2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleOffset(1.2);
    
    // determine draw type (TODO: probably could be generalized somehow)
    TString drawStr = "";
    switch(hist->GetDimension()) {
      case 1:
        drawStr = "EX0 P";
        if(histName=="Q_xsec") {
          hist->GetXaxis()->SetRangeUser(1,10);
          hist->GetYaxis()->SetRangeUser(1e-6,1);
          drawStr = "E P";
        };
        break;
      case 2:
        drawStr = "COLZ";
        break;
      case 3:
        drawStr = "BOX";
        break;
    };

    hist->Draw(drawStr);

    canv->SetGrid(1,1);
    canv->SetLogx(H->GetHistConfig(histName)->logx);
    canv->SetLogy(H->GetHistConfig(histName)->logy);
    canv->SetLogz(H->GetHistConfig(histName)->logz);
    canv->SetBottomMargin(0.15);
    canv->SetLeftMargin(0.15);
    canv->SetRightMargin(0.15);
    canv->Print(pngDir+"/"+canvN+".png");
    outfile->cd("/");
    canv->Write();
    outfile->cd("/");
  } else if (hist4 != nullptr) {
    TString canvN = "canv_"+histName+"_"+H->GetSetName();
    TCanvas *canv = new TCanvas(canvN,canvN, dimx, dimy);

    hist4->GetWaxis()->SetLabelSize(0.06);
    hist4->GetXaxis()->SetLabelSize(0.06);
    hist4->GetWaxis()->SetTitleSize(0.06);
    hist4->GetXaxis()->SetTitleSize(0.06);
    hist4->GetWaxis()->SetTitleOffset(1.2);
    
    //canv->SetGrid(1,1);
    canv->SetLogx(H->GetHist4Config(histName)->logx);
    canv->SetLogy(H->GetHist4Config(histName)->logy);
    canv->SetLogz(H->GetHist4Config(histName)->logz);
    canv->SetBottomMargin(0.15);
    canv->SetLeftMargin(0.15);
    canv->SetRightMargin(0.15);
    hist4->Draw();

    canv->Print(pngDir+"/"+canvN+".png");
    outfile->cd("/");
    canv->Write();
    outfile->cd("/");
  } else {
    cerr << "Couldn't find histogram " << histName << std::endl;
  }

  /* // deprecated, for combining single plots; TODO if needed, make a separate method
  TH1 *histClone = (TH1*) hist->Clone();
  histClone->SetLineColor  (nsum<nsumMax?summaryColor[nsum]:kBlack);
  histClone->SetMarkerColor(nsum<nsumMax?summaryColor[nsum]:kBlack);
  histClone->SetMarkerStyle(nsum<nsumMax?summaryStyle[nsum]:kFullCircle);

  if(nsum==0) {
    summaryCanv = new TCanvas(
        "summaryCanv_"+histName,
        "summaryCanv_"+histName,
        dimx,dimy
        );
    summaryCanv->SetGrid(1,1);
    summaryCanv->SetLogx(H->GetHistConfig(histName)->logx);
    summaryCanv->SetLogy(H->GetHistConfig(histName)->logy);
    summaryCanv->SetLogz(H->GetHistConfig(histName)->logz);
    summaryCanv->SetBottomMargin(0.15);
    summaryCanv->SetLeftMargin(0.15);
  };
  summaryCanv->cd();
  histClone->Draw(drawStr+(nsum>0?" SAME":""));
  nsum++;
  */
};

//=========================================================================
/* ALGORITHM: draw histograms from different bins in their respective bins
on axis of bin variables, e.g. Q2 vs x.
*/
// -- process one histArr
void PostProcessor::DrawInBins(
    TString outName,    
    std::vector<std::vector<Histos*>>& histArr,
    TString histName,
    TString var1name, int nvar1, double var1low, double var1high, bool var1log,
    TString var2name, int nvar2, double var2low, double var2high, bool var2log,
    bool intgrid1, bool intgrid2, // grid option for small plots
    bool renormalize // if true, normalize histograms
){
  std::vector<std::vector<std::vector<Histos*>>> histArrList;
  histArrList.push_back(histArr); // build list of histArrs with the one histArr
  this->DrawInBins( // call the main algo below, which uses histArrList
      outName,    
      histArrList,
      histName,
      var1name, nvar1, var1low, var1high, var1log,
      var2name, nvar2, var2low, var2high, var2log,
      intgrid1, intgrid2,
      renormalize
      );
};
// -- process list of histArrs
void PostProcessor::DrawInBins(
    TString outName,    
    std::vector<std::vector<std::vector<Histos*>>>& histArrList,
    TString histName,
    TString var1name, int nvar1, double var1low, double var1high, bool var1log,
    TString var2name, int nvar2, double var2low, double var2high, bool var2log,
    bool intgrid1, bool intgrid2, // grid option for small plots
    bool renormalize // if true, normalize histograms (useful for drawing 1D hists with "SAME")
){
  // default values set for nvar1==nvar2
  int canvx = 1400;
  int canvy = 1200;
  double botmargin = 0.2;
  double leftmargin = 0.2;
  double xaxisy = 0.04;
  double xaxisx1 = 0.08;
  double xaxisx2 = 0.97;
  double yaxisx = 0.04;
  double yaxisy1 = 0.085;
  double yaxisy2 = 0.97;
  if(nvar1 > nvar2){
    // different canvas sizing/axis position for unequal binning
    canvx = 2200;
    canvy = 1400;
    xaxisx1 = 0.075;
    xaxisx2 = 0.975;
    yaxisy1 = 0.08;
  };

  TString canvN = "canv_"+outName+"_"+histName;
  TCanvas *canv = new TCanvas(canvN,canvN, canvx, canvy);
  TPad *mainpad = new TPad("mainpad", "mainpad", 0.07, 0.07, 0.98, 0.98);

  TLegend *leg;
  if(legendLabels.size()>0) {
    double legX0 = 0.1;
    double legY0 = 0.95;
    double legWidth = 0.15;
    double legRowHeight = 0.05;
    leg = new TLegend(
        legX0,
        legY0,
        legX0 + legWidth,
        legY0 - legRowHeight*legendLabels.size()
        );
  };

  mainpad->SetFillStyle(4000);
  mainpad->Divide(nvar1,nvar2,0,0);
  mainpad->Draw();
  TLine * lDIRC = new TLine(6,-1,6,1);
  TLine * lDIRClow = new TLine(0.5,-1,0.5,1);
  TLine * lmRICH = new TLine(2,-1,2,-4);
  TLine * lDRICH = new TLine(2.5,1,2.5,4);
  lDIRC->SetLineColor(kRed);
  lDIRClow->SetLineColor(kRed);
  lmRICH->SetLineColor(kRed);
  lDRICH->SetLineColor(kRed);
  // TH1* histArray[nvar1][nvar2]; // TODO: re-enable
  int drawpid = 0;
  outfile->cd("/");
  canv->Write();

  // get histograms from Histos 2D vector
  Histos *H;
  for(int i = 0; i < nvar1; i++){
    for(int j = 0; j < nvar2; j++){
      int count = 0;
      int dims;
      for(auto histArr : histArrList) {
        H = histArr[i][j];
        TH1 *hist = H->Hist(histName);
        dims = hist->GetDimension();
        // histArray[i][j] = hist;
        hist->SetTitle("");
        //hist->GetXaxis()->SetTitle("");
        //hist->GetYaxis()->SetTitle("");
        //hist->GetXaxis()->SetLabelSize(0);
        //hist->GetYaxis()->SetLabelSize(0);

        hist->SetLineWidth(3);
        switch(count) {
          case 0:
            hist->SetLineColor(kBlack);
            hist->SetFillColor(kGray);
            break;
          case 1:
            hist->SetLineColor(kAzure);
            hist->SetLineStyle(1);
            break;
          case 2:
            hist->SetLineColor(kRed);
            hist->SetLineStyle(7);
            break;
          case 3:
            hist->SetLineColor(kGreen+2);
            hist->SetLineStyle(9);
            break;
        }

        if(legendLabels.size()>0 && i==0 && j==0) {
          try { leg->AddEntry(hist,legendLabels[count],"LF"); }
          catch(const std::out_of_range &e) {
            cerr << "ERROR: legendLabels not filled correctly" << endl;
          };
        };

        if(count==0) {
          mainpad->cd((nvar2-j-1)*nvar1 + i + 1);
          gPad->SetLogx(H->GetHistConfig(histName)->logx);
          gPad->SetLogy(H->GetHistConfig(histName)->logy);
          gPad->SetLogz(H->GetHistConfig(histName)->logz);
          gPad->SetGridy(intgrid2);
          gPad->SetGridx(intgrid1);
        }

        // renormalize
        if(renormalize) {
          // hist->Scale(1/hist->GetMaximum());
          // hist->Scale(1/hist->GetEntries());
          hist->Scale(1/hist->Integral());
        }

        TString drawStr = "";
        switch(dims) {
          case 1:
            drawStr = "HIST MIN0";
            break;
          case 2:
            drawStr = "COLZ";
            break;
          case 3:
            drawStr = "BOX";
            break;
        };
        if(count>0) drawStr += " SAME";

        if( hist->GetEntries() > 0 ) {	
          hist->Draw(drawStr);
          if(drawpid){
            lDIRClow->Draw();
            lDIRC->Draw();
            lmRICH->Draw();
            lDRICH->Draw();
          }
          hist->GetXaxis()->SetLabelSize(0.04);
          hist->GetYaxis()->SetLabelSize(0.04);
          hist->GetXaxis()->SetTitleSize(0.05);
          hist->GetYaxis()->SetTitleSize(0.05);
          hist->GetXaxis()->SetTitleOffset(0.9);
          hist->GetXaxis()->SetLabelOffset(0.0005);
          if(dims==1) {
            hist->GetYaxis()->SetLabelSize(0.00); // suppress y-axis labels (since each subplot has its own scale)
          }
        }
        count++;
      }; // end for(histArrList)

      // some formatting for 1D plots drawn with "SAME"
      if(dims==1) {
        UnzoomVertical(gPad,"",true,H->GetHistConfig(histName)->logy); // allow all "SAME" plots to be visible; 3rd arg forces min to be 0
      }
    };
  };
  canv->cd();

  TPad *newpad1 = new TPad("newpad1","full pad",0,0,1,1);
  TPad *newpad2 = new TPad("newpad2","full pad",0,0,1,1);
  newpad1->SetFillStyle(4000);
  newpad1->Draw();
  newpad2->SetFillStyle(4000);
  newpad2->Draw();

  TString xopt, yopt;
  if(var1log) xopt = "GS";
  else xopt = "S";
  if(var2log) yopt = "GS";
  else yopt = "S";

  TGaxis *xaxis = new TGaxis(xaxisx1,xaxisy,xaxisx2,xaxisy,var1low,var1high,510,xopt);
  TGaxis *yaxis = new TGaxis(yaxisx,yaxisy1,yaxisx,yaxisy2,var2low,var2high,510,yopt);
  xaxis->SetTitle(var1name);
  xaxis->SetName("xaxis");
  xaxis->SetTitleSize(0.02);
  xaxis->SetTextFont(40);
  xaxis->SetLabelSize(0.02);
  xaxis->SetTickSize(0.02);

  yaxis->SetTitle(var2name);
  yaxis->SetTitleSize(0.02);
  yaxis->SetName("yaxis");
  yaxis->SetTextFont(40);
  yaxis->SetLabelSize(0.02);
  yaxis->SetTickSize(0.02);

  newpad1->cd();
  yaxis->Draw();
  newpad2->cd();
  xaxis->Draw();

  if(legendLabels.size()>0) leg->Draw();

  //  canv->Write();
  canv->Print(pngDir+"/"+canvN+".png");
  canv->Print(pngDir+"/"+canvN+".pdf");
  outfile->cd("/");
  canv->Write();
  // for(int i = 0; i <nvar1; i++){
  //   for(int j = 0; j < nvar2; j++){
  //     histArray[i][j]->Write();
  //   }
  // }
};


//=========================================================================
/* ALGORITHM: draw ratios of histograms from different bins in their respective bins
on axis of bin variables, e.g. Q2 vs x.
*/
// -- process one histArr
void PostProcessor::DrawRatioInBins(
    TString outName,    
    std::vector<std::vector<Histos*>>& histArr,
    TString histName,
    TString var1name, int nvar1, double var1low, double var1high, bool var1log,
    TString var2name, int nvar2, double var2low, double var2high, bool var2log,
    bool intgrid1, bool intgrid2, // grid option for small plots
    bool renormalize // if true, normalize histograms
){
  std::vector<std::vector<std::vector<Histos*>>> histArrList;
  histArrList.push_back(histArr); // build list of histArrs with the one histArr
  this->DrawRatioInBins( // call the main algo below, which uses histArrList
      outName,    
      histArrList,
      histName,
      var1name, nvar1, var1low, var1high, var1log,
      var2name, nvar2, var2low, var2high, var2log,
      intgrid1, intgrid2,
      renormalize
      );
};
// -- process list of histArrs
void PostProcessor::DrawRatioInBins(
    TString outName,    
    std::vector<std::vector<std::vector<Histos*>>>& histArrList,
    TString histName,
    TString var1name, int nvar1, double var1low, double var1high, bool var1log,
    TString var2name, int nvar2, double var2low, double var2high, bool var2log,
    bool intgrid1, bool intgrid2, // grid option for small plots
    bool renormalize // if true, normalize histograms (useful for drawing 1D hists with "SAME")
){
  // default values set for nvar1==nvar2
  int canvx = 1400;
  int canvy = 1200;
  double botmargin = 0.2;
  double leftmargin = 0.2;
  double xaxisy = 0.04;
  double xaxisx1 = 0.08;
  double xaxisx2 = 0.97;
  double yaxisx = 0.04;
  double yaxisy1 = 0.085;
  double yaxisy2 = 0.97;
  if(nvar1 > nvar2){
    // different canvas sizing/axis position for unequal binning
    canvx = 2200;
    canvy = 1400;
    xaxisx1 = 0.075;
    xaxisx2 = 0.975;
    yaxisy1 = 0.08;
  };

  TString canvN = "canv_"+outName+"_"+histName;
  TCanvas *canv = new TCanvas(canvN,canvN, canvx, canvy);
  TPad *mainpad = new TPad("mainpad", "mainpad", 0.07, 0.07, 0.98, 0.98);

  TLegend *leg;
  if(legendLabels.size()>0) {
    double legX0 = 0.1;
    double legY0 = 0.95;
    double legWidth = 0.15;
    double legRowHeight = 0.05;
    leg = new TLegend(
        legX0,
        legY0,
        legX0 + legWidth,
        legY0 - legRowHeight*legendLabels.size()
        );
  };

  mainpad->SetFillStyle(4000);
  mainpad->Divide(nvar1,nvar2,0,0);
  mainpad->Draw();
  TLine * lDIRC = new TLine(6,-1,6,1);
  TLine * lDIRClow = new TLine(0.5,-1,0.5,1);
  TLine * lmRICH = new TLine(2,-1,2,-4);
  TLine * lDRICH = new TLine(2.5,1,2.5,4);
  lDIRC->SetLineColor(kRed);
  lDIRClow->SetLineColor(kRed);
  lmRICH->SetLineColor(kRed);
  lDRICH->SetLineColor(kRed);
  // TH1* histArray[nvar1][nvar2]; // TODO: re-enable
  int drawpid = 0;
  outfile->cd("/");
  canv->Write();

  // get histograms from Histos 2D vector
  Histos *H;
  for(int i = 0; i < nvar1; i++){
    for(int j = 0; j < nvar2; j++){
      int count = 0;
      int dims;
      TH1 *histref;
      for(auto histArr : histArrList) {
        H = histArr[i][j];
        TH1 *hist = H->Hist(histName);
        dims = hist->GetDimension();
        // histArray[i][j] = hist;
        hist->SetTitle("");
        //hist->GetXaxis()->SetTitle("");
        //hist->GetYaxis()->SetTitle("");
        //hist->GetXaxis()->SetLabelSize(0);
        //hist->GetYaxis()->SetLabelSize(0);

        hist->SetLineWidth(3);
        switch(count) {
          case 0:
	    histref = (TH1*) hist->Clone();
            hist->SetLineColor(kBlack);
	    //            hist->SetFillColor(kGray);
            break;
          case 1:
            hist->SetLineColor(kAzure);
            hist->SetLineStyle(1);
            break;
          case 2:
            hist->SetLineColor(kRed);
            hist->SetLineStyle(7);
            break;
          case 3:
            hist->SetLineColor(kGreen+2);
            hist->SetLineStyle(9);
            break;
        }


	hist->Divide(histref);
	
        if(legendLabels.size()>0 && i==0 && j==0) {
          try { leg->AddEntry(hist,legendLabels[count],"LF"); }
          catch(const std::out_of_range &e) {
            cerr << "ERROR: legendLabels not filled correctly" << endl;
          };
        };

        if(count==0) {
          mainpad->cd((nvar2-j-1)*nvar1 + i + 1);
          gPad->SetLogx(H->GetHistConfig(histName)->logx);
          gPad->SetLogy(H->GetHistConfig(histName)->logy);
          gPad->SetLogz(H->GetHistConfig(histName)->logz);
          gPad->SetGridy(intgrid2);
          gPad->SetGridx(intgrid1);
        }

        // renormalize
        if(renormalize) {
          // hist->Scale(1/hist->GetMaximum());
          // hist->Scale(1/hist->GetEntries());
	  //         hist->Scale(1/hist->Integral());
	  hist->GetYaxis()->SetRangeUser(0.,3.);
        }

        TString drawStr = "";
        switch(dims) {
          case 1:
            drawStr = "HIST MIN0";
            break;
          case 2:
            drawStr = "COLZ";
            break;
          case 3:
            drawStr = "BOX";
            break;
        };
        if(count>0) drawStr += " SAME";

        if( hist->GetEntries() > 0 ) {	
          hist->Draw(drawStr);
          if(drawpid){
            lDIRClow->Draw();
            lDIRC->Draw();
            lmRICH->Draw();
            lDRICH->Draw();
          }
          hist->GetXaxis()->SetLabelSize(0.04);
          hist->GetYaxis()->SetLabelSize(0.04);
          hist->GetXaxis()->SetTitleSize(0.05);
          hist->GetYaxis()->SetTitleSize(0.05);
          hist->GetXaxis()->SetTitleOffset(0.9);
          hist->GetXaxis()->SetLabelOffset(0.0005);
          if(dims==1) {
            hist->GetYaxis()->SetLabelSize(0.00); // suppress y-axis labels (since each subplot has its own scale)
          }
        }
        count++;
      }; // end for(histArrList)

      // some formatting for 1D plots drawn with "SAME"
      if(dims==1) {
        UnzoomVertical(gPad,"",true,H->GetHistConfig(histName)->logy); // allow all "SAME" plots to be visible; 3rd arg forces min to be 0
      }
    };
  };
  canv->cd();

  TPad *newpad1 = new TPad("newpad1","full pad",0,0,1,1);
  TPad *newpad2 = new TPad("newpad2","full pad",0,0,1,1);
  newpad1->SetFillStyle(4000);
  newpad1->Draw();
  newpad2->SetFillStyle(4000);
  newpad2->Draw();

  TString xopt, yopt;
  if(var1log) xopt = "GS";
  else xopt = "S";
  if(var2log) yopt = "GS";
  else yopt = "S";

  TGaxis *xaxis = new TGaxis(xaxisx1,xaxisy,xaxisx2,xaxisy,var1low,var1high,510,xopt);
  TGaxis *yaxis = new TGaxis(yaxisx,yaxisy1,yaxisx,yaxisy2,var2low,var2high,510,yopt);
  xaxis->SetTitle(var1name);
  xaxis->SetName("xaxis");
  xaxis->SetTitleSize(0.02);
  xaxis->SetTextFont(40);
  xaxis->SetLabelSize(0.02);
  xaxis->SetTickSize(0.02);

  yaxis->SetTitle(var2name);
  yaxis->SetTitleSize(0.02);
  yaxis->SetName("yaxis");
  yaxis->SetTextFont(40);
  yaxis->SetLabelSize(0.02);
  yaxis->SetTickSize(0.02);

  newpad1->cd();
  yaxis->Draw();
  newpad2->cd();
  xaxis->Draw();

  if(legendLabels.size()>0) leg->Draw();

  //  canv->Write();
  canv->Print(pngDir+"/"+canvN+".png");
  canv->Print(pngDir+"/"+canvN+".pdf");
  outfile->cd("/");
  canv->Write();
  // for(int i = 0; i <nvar1; i++){
  //   for(int j = 0; j < nvar2; j++){
  //     histArray[i][j]->Write();
  //   }
  // }
};


//=========================================================================

/* ALGORITHM: draw a ratio of all 1D histograms in the specified histogram set
* - the ratio will be of `numerSet` over `denomSet`
* - canvases will be created, and depending on the setting of `plotRatioOnly`,
*   either just the ratio will be plotted, or the two histograms will also be
*   plotted
* - Use `outName` to specify a names for output canvases
* - summary canvases are also created, which can combine multiple ratio plots;
*   they are accumulated in `summaryCanvMap` and colors are set by `nsum`; call
*   `ResetVars` to reset the accumulation
 * - when done looping, call `FinishDrawRatios`
*/
void PostProcessor::DrawRatios(
    TString outName, Histos *numerSet, Histos *denomSet, Bool_t plotRatioOnly
) {

  cout << "draw ratios " << outName << "..." << endl;
  enum HHenum {num,den};
  Histos *HH[2];
  HH[num] = numerSet;
  HH[den] = denomSet;
  TH1 *hist[2];
  TH1D *ratio;
  Double_t ny,ry,err;
  TCanvas *canv;
  TString histT[2], tok[2], uniqT[2], sameT[2];
  Ssiz_t tf;

  outfile->cd("/");
  outfile->mkdir(outName);
  outfile->cd(outName);

  // loop over 1D histograms
  for(TString varName : HH[num]->VarNameList) {
    hist[num] = HH[num]->Hist(varName);
    if(hist[num] != nullptr && hist[num]->GetDimension()==1) {
      hist[den] = HH[den]->Hist(varName);

      // filter title
      // TODO: there is a bug here now, ratio denominators aren't quite right
      for(int f=0; f<2; f++)  histT[f] = hist[f]->GetTitle();
      for(int f=0; f<2; f++) {
        tf = 0;
        sameT[f] = "";
        uniqT[f] = "";
        while(histT[f].Tokenize(tok[f],tf,", ")) {
          if(histT[(f+1)%2].Contains(tok[f])) sameT[f] += ", " + tok[f];
          else                                uniqT[f] += ", " + tok[f];
        };
        sameT[f](TRegexp("^, "))="";
        uniqT[f](TRegexp("^, "))="";
        /*cout <<  "histT = " << histT[f] << endl;
          cout <<  "uniqT = " << uniqT[f] << endl;
          cout <<  "sameT = " << sameT[f] << endl;*/
        hist[f]->SetTitle(sameT[f]);
      };
      TString ratioStr = "ratio("+uniqT[num]+"/"+uniqT[den]+")";

      // calculate ratio
      ratio = (TH1D*) hist[num]->Clone();
      ratio->Divide(hist[den]);
      TString ratioT = ratio->GetTitle();
      TString ratioSummaryT = ratioT;
      ratioT(TRegexp("distribution")) = ratioStr;
      ratioSummaryT(TRegexp("distribution")) = "ratios";
      ratio->SetTitle(ratioT);

      // calculate ratio errors, assuming ratio represents an efficiency
      for(int b=0; b<=ratio->GetNbinsX(); b++) {
        ny = hist[den]->GetBinContent(b);
        ry = ratio->GetBinContent(b);
        err = TMath::Abs(ny)>0 ?
          TMath::Sqrt( ry*(1-ry) / ny )
          : 0;
        ratio->SetBinError(b,err);
      };

      // canvas
      TString canvN = "canv_"+outName+"_"+varName;
      Int_t dimx=800;
      Int_t dimy=700;
      canv = new TCanvas(canvN,canvN, dimx*(plotRatioOnly?1:2), dimy);

      // draw plots
      int numPads;
      if(!plotRatioOnly) { canv->Divide(2,1); numPads=2; }
      else { canv->Divide(1,1); numPads=1; };
      if(!plotRatioOnly) {
        canv->cd(2);
        hist[num]->SetLineColor(kYellow-8);
        hist[den]->SetLineColor(kAzure-7);
        hist[num]->SetMarkerColor(kYellow-8);
        hist[den]->SetMarkerColor(kAzure-7);
        hist[num]->SetMarkerStyle(kFullCircle);
        hist[den]->SetMarkerStyle(kFullCircle);
        for(int f=0; f<2; f++) {
          hist[f]->SetMarkerSize(1.0);
          hist[f]->SetLineWidth(2);
          hist[f]->GetXaxis()->SetLabelSize(0.06);
          hist[f]->GetYaxis()->SetLabelSize(0.06);
          hist[f]->GetXaxis()->SetTitleSize(0.06);
          hist[f]->GetYaxis()->SetTitleSize(0.06);
          hist[f]->GetXaxis()->SetTitleOffset(1.2);
          //if(plotName=="PperpDistLin") hist[f]->GetXaxis()->SetRangeUser(0,1.7);
        };
        hist[den]->Draw("EX0 P");
        hist[num]->Draw("EX0 P SAME");
      };
      canv->cd(1);
      ratio->GetYaxis()->SetTitle(ratioStr);
      ratio->SetLineColor(kBlack);
      ratio->SetMarkerColor(kBlack);
      ratio->SetMarkerStyle(kFullCircle);
      ratio->SetMarkerSize(1.0);
      ratio->SetLineWidth(2);
      ratio->GetXaxis()->SetLabelSize(0.06);
      ratio->GetYaxis()->SetLabelSize(0.06);
      ratio->GetXaxis()->SetTitleSize(0.06);
      ratio->GetYaxis()->SetTitleSize(0.06);
      ratio->GetXaxis()->SetTitleOffset(1.2);
      ratio->GetYaxis()->SetRangeUser(0,1.3);
      ratio->Draw("EX0 P");
      for(int pad=1; pad<=numPads; pad++) {
        canv->GetPad(pad)->SetGrid(1,1);
        canv->GetPad(pad)->SetLogx(HH[num]->GetHistConfig(varName)->logx);
        if(pad>1) {
          canv->GetPad(pad)->SetLogy(HH[num]->GetHistConfig(varName)->logy);
          canv->GetPad(pad)->SetLogz(HH[num]->GetHistConfig(varName)->logz);
        };
        canv->GetPad(pad)->SetBottomMargin(0.15);
        canv->GetPad(pad)->SetLeftMargin(0.15);
      };
      canv->Print(pngDir+"/"+canvN+".png");
      canv->Write();

      // add to summary canvas
      TH1D *ratioSummary = (TH1D*) ratio->Clone();
      ratioSummary->SetTitle(ratioSummaryT);
      ratioSummary->SetYTitle("ratio");
      ratioSummary->SetLineColor  (nsum<nsumMax?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerColor(nsum<nsumMax?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerStyle(nsum<nsumMax?summaryStyle[nsum]:kFullCircle);
      auto kv = summaryCanvMap.find(varName);
      if(kv==summaryCanvMap.end()) {
        summaryCanv = new TCanvas(
            "summaryCanv_ratio"+varName,
            "summaryCanv_ratio"+varName,
            dimx,dimy
            );
        summaryCanv->SetGrid(1,1);
        summaryCanv->SetLogx(HH[num]->GetHistConfig(varName)->logx);
        summaryCanv->SetBottomMargin(0.15);
        summaryCanv->SetLeftMargin(0.15);
        summaryCanvMap.insert({varName,summaryCanv});
        ratioSummary->Draw();
      } else {
        kv->second->cd();
        ratioSummary->Draw("SAME");
      };

    };
  };
  outfile->cd("/");
  nsum++;
};
void PostProcessor::FinishDrawRatios(TString summaryDir) {
  std::cout << "CALL PostProcessor::FinishDrawRatios" << endl;
  // write summary canvas
  outfile->cd("/");
  outfile->mkdir(summaryDir);
  outfile->cd(summaryDir);
  for(auto const& kv : summaryCanvMap) {
    kv.second->Write();
    kv.second->Print(pngDir+"/"+summaryDir+"_"+kv.first+".png");
  };
  outfile->cd("/");
  this->ResetVars();
};



//=========================================================================
// TEXT FILE: start new text file, append to text file, print text file
void PostProcessor::StartTextFile(TString datFile, TString firstLine) {
  gSystem->RedirectOutput(datFile,"w");
  if(firstLine!="") cout << firstLine << endl;
  gSystem->RedirectOutput(0);
};
void PostProcessor::AppendToTextFile(TString datFile, TString appendText="") {
  gSystem->RedirectOutput(datFile,"a");
  cout << appendText << endl;
  gSystem->RedirectOutput(0);
};
void PostProcessor::PrintTextFile(TString datFile) {
  gROOT->ProcessLine(".! cat "+datFile);
};


//=========================================================================
// TEXT FILE: pipe `inputFile` through `column -t` and append output to `outputFile`
void PostProcessor::Columnify(TString inputFile, TString outputFile) {
  gSystem->RedirectOutput(outputFile,"a");
  gROOT->ProcessLine(".! cat "+inputFile+" | column -t");
  gSystem->RedirectOutput(0);
};


//=========================================================================
// get pointer to BinSet associated with variable
BinSet *PostProcessor::GetBinSet(TString varName) {
  return (BinSet*) infile->Get(varName+"_bins");
};
// get pointer to CutDef associated with a particular bin
CutDef *PostProcessor::GetBinCut(TString varName, Int_t binNum) {
  return this->GetBinSet(varName)->Cut(binNum);
};
// get list of bin numbers
std::vector<int> PostProcessor::GetBinNums(TString varName) {
  BinSet *B = this->GetBinSet(varName);
  std::vector<int> retVec;
  if(B==nullptr) {
    cerr << "ERROR: unknown variable " << varName << " in PostProcessor::GetBinNums" << endl;
  } else {
    for(int n=0; n<B->GetNumBins(); n++) retVec.push_back(n);
  };
  return retVec;
};

// return true if the bin is "full" range, and it's not the only bin
Bool_t PostProcessor::SkipFull(TString varName, Int_t binNum) {
  return GetBinSet(varName)->GetNumBins() > 1 &&
         GetBinCut(varName,binNum)->GetCutType()=="Full";
};


//=========================================================================
// destructor
PostProcessor::~PostProcessor() {
};

