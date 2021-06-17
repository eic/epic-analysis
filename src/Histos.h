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

    // builders
    void DefineHist(
        TString varname,
        TString vartitle,
        TString units,
        Int_t nb, Double_t lb, Double_t ub
        );

    // accessors
    TH1 *Hist(TString histName);

    // writers
    void WriteHists() {
      histList->Write(setname+"_hists",TObject::kSingleKey);
    };


  private:
    TString setname,settitle;
    TObjArray *histList;
    map<TString,TH1*> histMap;

  ClassDef(Histos,1);
};

#endif
