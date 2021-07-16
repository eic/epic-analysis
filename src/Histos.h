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
#include "TMath.h"

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

class Histos : public TNamed
{
  public:
    Histos(TString setname_="setname", TString settitle_="settitle");
    ~Histos();

    // number of bins
    const Int_t NBINS = 50;

    // accessors
    TH1 *Hist(TString histName); // access histogram by name
    vector<TString> VarNameList; // list of histogram names (for external looping)

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

    // writers
    void WriteHists() {
      histList->Write("histArr_"+setname,TObject::kSingleKey);
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
