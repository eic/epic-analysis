// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Duane Byer

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

// ROOT
#include <TSystem.h>
#include <TObject.h>
#include <TNamed.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>

// adage
#include <adage/BinSet.h>
#include <adage/CutDef.h>

// epic-analysis
#include "Hist4D.h"

// container for histogram settings
class HistConfig : public TNamed {
  public:
    Bool_t logx;
    Bool_t logy;
    Bool_t logz;
    Bool_t logw;
    HistConfig() {
      logx=false;
      logy=false;
      logz=false;
      logw=false;
    };
    ~HistConfig() {};
  ClassDef(HistConfig,1);
};

// container for histograms
class Histos : public TNamed
{
  public:
    Histos(TString setname_="setname", TString settitle_="settitle");
    ~Histos();

    // accessors
    TH1 *Hist(TString histName, Bool_t silence=false); // access histogram by name
    Hist4D *Hist4(TString histName, Bool_t silence=false);
    HistConfig *GetHistConfig(TString histName); // settings for this histogram
    HistConfig *GetHist4Config(TString histName);
    std::vector<TString> VarNameList; // list of histogram names (for external looping)
    std::vector<CutDef*> CutDefList; // list of associated cut definitions
    TString GetSetName() { return setname; };
    TString GetSetTitle() { return settitle; };
    CutDef *GetCutDef(TString varName); // get a specific CutDef

    // store associated cut definitions
    void AddCutDef(CutDef *cut) { CutDefList.push_back(cut); };

    // histogram builders
    void DefineHist1D(
        TString varname,
        TString vartitle,
        TString units,
        Int_t numBins, Double_t lowerBound, Double_t upperBound,
        Bool_t logx = false,
        Bool_t logy = false
        );
    void DefineHist2D(
        TString varname,
        TString vartitlex, TString vartitley,
        TString unitsx, TString unitsy,
        Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
        Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
        Bool_t logx = false,
        Bool_t logy = false,
        Bool_t logz = false
        );
    void DefineHist2D(
        TString varname,
	TString vartitlex, TString vartitley,
	TString unitsx, TString unitsy,
	Int_t numBinsx, Double_t *xBins,
	Int_t numBinsy, Double_t *yBins,
	Bool_t logx = false,
	Bool_t logy = false,
	Bool_t logz = false
	);

    void DefineHist3D(
        TString varname,
        TString vartitlex, TString vartitley, TString vartitlez,
        TString unitsx, TString unitsy, TString unitsz,
        Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
        Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
        Int_t numBinsz, Double_t lowerBoundz, Double_t upperBoundz,
        Bool_t logx = false,
        Bool_t logy = false,
        Bool_t logz = false
        );
    void DefineHist4D(
        TString varname,
        TString vartitlew, TString vartitlex, TString vartitley, TString vartitlez,
        TString unitsw, TString unitsx, TString unitsy, TString unitsz,
        Int_t numBinsw, Double_t lowerBoundw, Double_t upperBoundw,
        Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
        Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
        Int_t numBinsz, Double_t lowerBoundz, Double_t upperBoundz,
        Bool_t logw = false,
        Bool_t logx = false,
        Bool_t logy = false,
        Bool_t logz = false
        );

    // writers
    void WriteHists(TFile *ofile) {
      ofile->cd("/");
      ofile->mkdir("histArr_"+setname);
      ofile->cd("histArr_"+setname);
      for(auto const &kv : histMap) kv.second->Write();
      for(auto const &kv : hist4Map) kv.second->Write();
      ofile->cd("/");
    };

    // add histograms of `in_histos` to histograms in `this`
    void AddHistos(Histos *in_histos);

    // fill a histogram, iff it is defined
    template <typename... VALS> void FillHist1D(TString histName, VALS... values) {
      auto hist = Hist(histName,true);
      if(hist) hist->Fill(values...);
    }
    template <typename... VALS> void FillHist2D(TString histName, VALS... values) {
      auto hist = Hist(histName,true);
      if(hist) dynamic_cast<TH2*>(hist)->Fill(values...);
    }
    template <typename... VALS> void FillHist3D(TString histName, VALS... values) {
      auto hist = Hist(histName,true);
      if(hist) dynamic_cast<TH3*>(hist)->Fill(values...);
    }
    template <typename... VALS> void FillHist4D(TString histName, VALS... values) {
      auto hist = Hist4(histName,true);
      if(hist) hist->Fill(values...);
    }

  private:
    TString setname,settitle;
    std::map<TString,TH1*> histMap;
    std::map<TString,Hist4D*> hist4Map;
    std::map<TString,HistConfig*> histConfigMap;
    std::map<TString,HistConfig*> hist4ConfigMap;
    void RegisterHist(TString varname_, TH1 *hist_, HistConfig *config_);
    void RegisterHist4(TString varname_, Hist4D *hist_, HistConfig *config_);

  ClassDef(Histos,1);
};
