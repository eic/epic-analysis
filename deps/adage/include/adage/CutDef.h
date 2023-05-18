#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include <TSystem.h>
#include <TObject.h>
#include <TNamed.h>
#include <TString.h>
#include <TMath.h>


class CutDef : public TObject
{
  public:
    CutDef();
    // primary constructor (see CutDef.cxx)
    CutDef(
        TString varName_, TString varTitle_, TString cutType_,
        Double_t arg1=-1, Double_t arg2=-1
        );
    // constructor for just storing a cut ID (for externally applied cuts)
    CutDef(
        TString varName_, TString varTitle_,
        TString cutID_, TString cutTitle_
        );
    ~CutDef();

    // apply cut
    Bool_t CheckCut(Double_t arg1=-1);

    // accessors
    TString GetCutTitle() { return cutTitle; };
    TString GetVarName() { return varName; };
    TString GetVarTitle() { return varTitle; };
    TString GetCutType() { return cutType; };
    Double_t GetMin() { return min; };
    Double_t GetMax() { return max; };
    TString GetCutID() { return cutID; };
    Bool_t IsExternal() { return cutType.CompareTo("External",TString::kIgnoreCase)==0; };


  private:
    TString varName,varTitle,cutType;
    TString cutTitle;
    TString cutID;
    Double_t min,max;
    Double_t center,delta;

  ClassDef(CutDef,1);
};
