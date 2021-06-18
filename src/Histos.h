#ifndef Histos_
#define Histos_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using std::map;
using std::cout;
using std::cerr;
using std::endl;

class Histos : public TObject
{
  public:
    Histos(TString setname_, TString settitle_);
    ~Histos();

    const Int_t NBINS = 100;

    // histogram builders
    void DefineHist1D(
        TString varname,
        TString vartitle,
        TString units,
        Int_t numBins, Double_t lowerBound, Double_t upperBound,
        Bool_t logx = false
        );
    void DefineHist2D(
        TString varname,
        TString vartitlex, TString vartitley,
        TString unitsx, TString unitsy,
        Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
        Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
        Bool_t logx = false,
        Bool_t logy = false
        );

    // accessors
    TH1 *Hist(TString histName);

    // writers
    void WriteHists() {
      histList->Write(setname+"_hists",TObject::kSingleKey);
    };

    // tools
    static void BinLog(TAxis *axis);


  private:
    TString setname,settitle;
    TObjArray *histList;
    map<TString,TH1*> histMap;

  ClassDef(Histos,1);
};

#endif
