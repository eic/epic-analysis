#include <map>

// root
#include "TFile.h"
#include "TCanvas.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"

// largex-eic
#include "Histos.h"
#include "Kinematics.h"

// globals
TString infileN, outfileN, pngDir;
TFile *infile, *outfile;
Bool_t plotRatioOnly;
std::map<TString,TCanvas*> summaryCanvMap;
TCanvas *summaryCanv;
int nsum;


// subroutines
void DrawRatios(TString outName, TString numerSet, TString denomSet);

//=========================================================================

int main(int argc, char **argv) {

  // ARGUMENTS ////////////////////////////////////////////////
  infileN="out/histos.root";
  plotRatioOnly=false;
  if(argc>1) infileN = TString(argv[1]);
  if(argc>2) plotRatioOnly = ((Int_t)strtof(argv[2],NULL))>0;
  /////////////////////////////////////////////////////////////

  cout << "-- drawing histograms from " << infileN << endl;
  outfileN = infileN;
  outfileN(TRegexp("\\.root$")) = "";
  pngDir = outfileN+".images";
  outfileN += ".canvas.root";
  outfile = new TFile(outfileN,"RECREATE");
  infile = new TFile(infileN,"READ");

  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".! mkdir -p "+pngDir);
  nsum=0;

  // ratios y>y_min / y>0
  for(int y=1; y<4; y++) {
    DrawRatios(
        Form("yrat%d",y),
        Form("histos_pipTrack_ycut%d",y),
        "histos_pipTrack_ycut0"
        );
  };

  // write summary canvases
  outfile->cd("/");
  outfile->mkdir("summary");
  outfile->cd("summary");
  for(auto const& kv : summaryCanvMap) {
    kv.second->Write();
    kv.second->Print(pngDir+"/summary_"+kv.first+".png");
  };


  // cleanup
  infile->Close();
  outfile->Close();
  cout << outfileN << " written." << endl;
  cout << pngDir << "/ images written." << endl;

};


//=========================================================================

void DrawRatios(TString outName, TString numerSet, TString denomSet) {

  cout << "draw ratios " << outName << "..." << endl;
  enum HHenum {num,den};
  Histos *HH[2];
  HH[num] = (Histos*) infile->Get(numerSet);
  HH[den] = (Histos*) infile->Get(denomSet);
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
    if(hist[num]->GetDimension()==1) {
      hist[den] = HH[den]->Hist(varName);

      // filter title
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
        canv->cd(1);
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
        canv->cd(2);
      } else {
        canv->cd(1);
      };
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
        canv->GetPad(pad)->SetLogy(HH[num]->GetHistConfig(varName)->logy);
        canv->GetPad(pad)->SetLogz(HH[num]->GetHistConfig(varName)->logz);
        canv->GetPad(pad)->SetBottomMargin(0.15);
        canv->GetPad(pad)->SetLeftMargin(0.15);
      };
      canv->Print(pngDir+"/"+canvN+".png");
      canv->Write();

      // add to summary canvas
      TH1D *ratioSummary = (TH1D*) ratio->Clone();
      ratioSummary->SetTitle(ratioSummaryT);
      ratioSummary->SetYTitle("ratio");
      Color_t summaryColor[3] = {kRed,kGreen+1,kBlue};
      Style_t summaryStyle[3] = {kFullCircle,kFullTriangleUp,kFullTriangleDown};
      ratioSummary->SetLineColor  (nsum<3?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerColor(nsum<3?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerStyle(nsum<3?summaryStyle[nsum]:kFullCircle);
      auto kv = summaryCanvMap.find(varName);
      if(kv==summaryCanvMap.end()) {
        summaryCanv = new TCanvas(
            "summaryCanv_"+varName,
            "summaryCanv_"+varName,
            dimx,dimy
            );
        summaryCanv->SetGrid(1,1);
        summaryCanv->SetLogx(HH[num]->GetHistConfig(varName)->logx);
        summaryCanv->SetBottomMargin(0.15);
        summaryCanv->SetLeftMargin(0.15);
        summaryCanvMap.insert(std::pair<TString,TCanvas*>(varName,summaryCanv));
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
