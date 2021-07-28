#ifndef Histos_
#define Histos_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"

// largex-eic
#include "BinSet.h"
#include "CutDef.h"

// container for histogram settings
class HistConfig : public TNamed {
  public:
    Bool_t logx;
    Bool_t logy;
    Bool_t logz;
    HistConfig() {
      logx=false;
      logy=false;
      logz=false;
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
    TH1 *Hist(TString histName); // access histogram by name
    HistConfig *GetHistConfig(TString histName); // settings for this histogram
    std::vector<TString> VarNameList; // list of histogram names (for external looping)
    std::vector<CutDef*> CutDefList; // list of associated cut definitions
    TString GetSetName() { return setname; };
    TString GetSetTitle() { return settitle; };

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

    // writers
    void WriteHists(TFile *ofile) {
      ofile->cd("/");
      ofile->mkdir("histArr_"+setname);
      ofile->cd("histArr_"+setname);
      for(auto const &kv : histMap) kv.second->Write();
      ofile->cd("/");
    };


  private:
    TString setname,settitle;
    std::map<TString,TH1*> histMap;
    std::map<TString,HistConfig*> histConfigMap;
    void RegisterHist(TString varname_, TH1 *hist_, HistConfig *config_);

  ClassDef(Histos,1);
};

#endif
