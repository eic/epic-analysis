#include "BinSet.h"

ClassImp(BinSet)

using std::cout;
using std::cerr;
using std::endl;

BinSet::BinSet(
    TString varName_, TString varTitle_,
    Int_t nbins_, Double_t min_, Double_t max_, Bool_t log_
    )
  : varName(varName_)
  , varTitle(varTitle_)
  , nbins(nbins_)
  , min(min_)
  , max(max_)
  , log(log_)
{
  // define axis
  axis = new TAxis(nbins,min,max);

  // logarithmic binning
  if(log) BinLog(axis);

  // build `bins` container
  bins = new TObjArray();
  for(int b=1; b<=axis->GetNbins(); b++) {
    bins->AddLast(
        new CutDef(
          varName,
          varTitle,
          "Range",
          axis->GetBinLowEdge(b),
          axis->GetBinUpEdge(b) 
          )
        );
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
  Float_t * newBins = new Float_t[nb+1];
  for (int b=0; b<=nb; b++) newBins[b] = TMath::Power(10,lb+b*wd);
  ax->Set(nb,newBins);
  delete[] newBins;
}; 


BinSet::~BinSet() {
};

