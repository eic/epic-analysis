#include "CutDef.h"

ClassImp(CutDef)

using std::cout;
using std::cerr;
using std::endl;

// default constructor (for streaming compatibility only)
CutDef::CutDef()
  : varName("unknown")
  , varTitle("unknown")
  , cutType("Full")
  , min(-1)
  , max(-1)
  , center(-1)
  , delta(-1)
{};

// constructor (for usage)
CutDef::CutDef(
    TString varName_, TString varTitle_, TString cutType_,
    Double_t arg1, Double_t arg2
    ) 
  : varName(varName_)
  , varTitle(varTitle_)
  , cutType(cutType_)
  , min(-1)
  , max(-1)
  , center(-1)
  , delta(-1)
{

  // minimum cut
  if(cutType.CompareTo("Min",TString::kIgnoreCase)==0) {
    min = arg1;
    cutTitle = Form("%s>%.2f",varTitle.Data(),min);
  }

  // maximum cut
  else if(cutType.CompareTo("Max",TString::kIgnoreCase)==0) {
    max = arg1;
    cutTitle = Form("%s<%.2f",varTitle.Data(),max);
  }

  // range (arg1,arg2)
  else if(cutType.CompareTo("Range",TString::kIgnoreCase)==0) {
    min = arg1;
    max = arg2;
    cutTitle = Form("%.2f<%s<%.2f",min,varTitle.Data(),max);
  }

  // arg1 +/- arg2
  else if(cutType.CompareTo("CenterDelta",TString::kIgnoreCase)==0) {
    center = arg1;
    delta = arg2;
    cutTitle = Form("%s#in%.2f#pm%.2f",varTitle.Data(),center,delta);
  }

  // full (no) cut
  else if(cutType.CompareTo("Full",TString::kIgnoreCase)==0) {
    cutTitle = Form("full %s",varTitle.Data());
  }

  // default to no cut
  else {
    cerr << "ERROR: unknown cutType " << cutType << endl;
    cerr << "       default to no cut" << endl;
    cutTitle="ERROR";
    cutType="Full";
  };
  

};


// apply cut
Bool_t CutDef::Cut(Double_t arg1) {

  // minimum cut
  if(cutType.CompareTo("Min",TString::kIgnoreCase)==0) {
    return arg1 > min;
  }

  // maximum cut
  else if(cutType.CompareTo("Max",TString::kIgnoreCase)==0) {
    return arg1 < max;
  }

  // range (arg1,arg2)
  else if(cutType.CompareTo("Range",TString::kIgnoreCase)==0) {
    return arg1 > min &&
           arg1 < max;
  }

  // arg1 +/- arg2
  else if(cutType.CompareTo("CenterDelta",TString::kIgnoreCase)==0) {
    return TMath::Abs(arg1-center) < delta;
  }

  // full (no) cut
  else if(cutType.CompareTo("Full",TString::kIgnoreCase)==0) {
    return true;
  };

  return false;
};

CutDef::~CutDef() {
};

