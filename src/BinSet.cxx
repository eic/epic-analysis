#include "BinSet.h"

ClassImp(BinSet)

using std::cout;
using std::cerr;
using std::endl;

// constructor
BinSet::BinSet(TString varName_, TString varTitle_)
  : varName(varName_)
  , varTitle(varTitle_)
  , binList(new TObjArray())
{
};

// build single bin
void BinSet::BuildBin(TString cutType_, Double_t arg1_, Double_t arg2_) {
  return BuildBin(
      new CutDef(
        varName,
        varTitle,
        cutType_,
        arg1_,
        arg2_
        ));
};
void BinSet::BuildBin(CutDef *cut_) { binList->AddLast(cut_); };

// build list of bins
void BinSet::BuildBins(Int_t nbins_, Double_t min_, Double_t max_, Bool_t log_) {
  BuildBins(new TAxis(nbins_,min_,max_),log_);
};
void BinSet::BuildBins(TAxis *ax, Bool_t log_) {
  if(log_) BinLog(ax);
  for(int b=1; b<=ax->GetNbins(); b++) {
    binList->AddLast(
        new CutDef(
          varName,
          varTitle,
          "Range",
          ax->GetBinLowEdge(b),
          ax->GetBinUpEdge(b) 
          ));
  };
};


// make equal-width log-scale bins
void BinSet::BinLog(TAxis *ax) {
  Float_t lb = ax->GetXmin();
  Float_t ub = ax->GetXmax();
  if(lb<=0||ub<=0||lb>=ub) {
    fprintf(stderr,"ERROR: bad axis range for BinSet::BinLog\n");
    return;
  };
  lb = TMath::Log10(lb);
  ub = TMath::Log10(ub);
  Int_t nb = ax->GetNbins();
  Float_t wd = (ub-lb)/nb;
  Float_t *newBins = new Float_t[nb+1];
  for(int b=0; b<=nb; b++) newBins[b] = TMath::Power(10,lb+b*wd);
  ax->Set(nb,newBins);
  delete[] newBins;
}; 


BinSet::~BinSet() {
};

