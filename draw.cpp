// root
#include "TFile.h"
#include "TCanvas.h"
#include "TRegexp.h"

// largex-eic
#include "Histos.h"
#include "Kinematics.h"

// globals
TString infileN;
TFile *infile, *outfile;
Bool_t plotRatioOnly;

// subroutines
void DrawRatios(Histos *numerSet, Histos *denomSet);

//=========================================================================

int main(int argc, char **argv) {

  // ARGUMENTS ////////////////////////////////////////////////
  infileN="out/histos.root";
  plotRatioOnly=false;
  if(argc>1) infileN = TString(argv[1]);
  if(argc>2) plotRatioOnly = ((Int_t)strtof(argv[2],NULL))>0;
  /////////////////////////////////////////////////////////////

  TString outfileN = infileN;
  outfileN(TRegexp("\\.root$")) = ".canvas.root";
  outfile = new TFile(outfileN,"RECREATE");
  infile = new TFile(infileN,"READ");

  Histos *h= (Histos*)infile->Get("histos_pipTrack_ycut3");
  DrawRatios(
      (Histos*)infile->Get("histos_pipTrack_ycut3"),
      (Histos*)infile->Get("histos_pipTrack_ycut0")
      );

  infile->Close();
  outfile->Close();

};


//=========================================================================

void DrawRatios(Histos *numerSet, Histos *denomSet) {
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

  // loop over 1D histograms
  for(TString varname : HH[num]->VarNameList) {
    hist[num] = HH[num]->Hist(varname);
    if(hist[num]->GetDimension()==1) {
      hist[den] = HH[den]->Hist(varname);

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
        cout <<  "histT = " << histT[f] << endl;
        cout <<  "uniqT = " << uniqT[f] << endl;
        cout <<  "sameT = " << sameT[f] << endl;
        hist[f]->SetTitle(sameT[f]);
      };
      TString ratioStr = "ratio("+uniqT[num]+"/"+uniqT[den]+")";

      // calculate ratio
      ratio = (TH1D*) hist[num]->Clone();
      ratio->Divide(hist[den]);
      TString ratioT = ratio->GetTitle();
      ratioT(TRegexp("distribution")) = ratioStr;
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
      canv = new TCanvas("canv_"+varname,"canv_"+varname, plotRatioOnly ? 800 : 1600, 700);

      if(!plotRatioOnly) {
        canv->Divide(2,1);
        canv->cd(1);
        canv->GetPad(1)->SetGrid(1,1);
        //canv->GetPad(1)->SetLogx();
        //canv->GetPad(1)->SetLogy();
        canv->GetPad(1)->SetBottomMargin(0.15);
        canv->GetPad(1)->SetLeftMargin(0.15);
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
        hist[den]->Draw("ERR P");
        hist[num]->Draw("ERR P SAME");
        canv->cd(2);
        canv->GetPad(2)->SetGrid(1,1);
        canv->GetPad(2)->SetBottomMargin(0.15);
        canv->GetPad(2)->SetLeftMargin(0.15);
      } else {
        canv->cd();
        canv->SetGrid(1,1);
        canv->SetBottomMargin(0.15);
        canv->SetLeftMargin(0.15);
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
      //if(plotName=="PperpDistLin") ratio->GetXaxis()->SetRangeUser(0,1.7);
      ratio->Draw("ERR X0 P");
      //canv->Print(Form("%s_bin%d.png",outfileN.Data(),b),"png");
      outfile->cd();
      canv->Write();
    };
  };
};
