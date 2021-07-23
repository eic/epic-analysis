#ifndef CutDef_
#define CutDef_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TMath.h"


class CutDef : public TObject
{
  public:
    CutDef();
    CutDef(
        TString varName_, TString varTitle_, TString cutType_,
        Double_t arg1=-1, Double_t arg2=-1
        );
    ~CutDef();

    // apply cut
    Bool_t Cut(Double_t arg1=-1);

    // accessors
    TString GetCutTitle() { return cutTitle; };
    TString GetVarName() { return varName; };


  private:
    TString varName,varTitle,cutType;
    TString cutTitle;
    Double_t min,max;
    Double_t center,delta;

  ClassDef(CutDef,1);
};

#endif
