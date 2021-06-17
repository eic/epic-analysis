#include "Histos.h"

ClassImp(Histos)


Histos::Histos(TString setname_, TString settitle_) {
  setname=setname_;
  settitle=settitle_;
  histList = new TObjArray();
  DefineHist("p","p","GeV",NBINS,0,10);
  DefineHist("pT","p_{T}^{lab}","GeV",NBINS,0,10);
  DefineHist("eta","#eta","",NBINS,-5,5);
  DefineHist("phi","#phi","",NBINS,-TMath::Pi(),TMath::Pi());
};


void Histos::DefineHist(
    TString varname,
    TString vartitle,
    TString units,
    Int_t nb, Double_t lb, Double_t ub
    ) {
  if(units!="") units=" ["+units+"]";
  TH1D *hist = new TH1D(
      setname+"_hist_"+varname,
      settitle+" "+vartitle+" distribution;"+vartitle+units,
      nb,lb,ub
      );
  histList->AddLast(hist);
  histMap.insert(std::pair<TString,TH1D*>(varname,hist));
};


TH1 *Histos::Hist(TString histName) {
  TH1 *retHist;
  try { retHist = histMap.at(histName); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: histMap does not have " 
         << histName << "histogram" << endl;
    return nullptr;
  };
  return retHist;
};


Histos::~Histos() {
};

