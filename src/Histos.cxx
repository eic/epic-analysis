#include "Histos.h"

ClassImp(Histos)


Histos::Histos(TString setname_, TString settitle_) {
  setname=setname_;
  settitle=settitle_;
  this->SetName("histos_"+setname);
  histList = new TObjArray();

  // HISTOGRAMS ================================================
  // DIS kinematics
  DefineHist2D("Q2vsX","x","Q^{2}","","GeV^{2}",
      NBINS,1e-3,1,
      NBINS,1e-3,100,
      true,true
      );
  DefineHist1D("W","W","GeV",NBINS,0,15);
  DefineHist1D("y","y","",NBINS,1e-5,1,true);
  // hadron 4-momentum
  DefineHist1D("p","p","GeV",NBINS,0,10);
  DefineHist1D("pTlab","p_{T}^{lab}","GeV",NBINS,1e-4,10,true);
  DefineHist1D("eta","#eta","",NBINS,-5,5);
  DefineHist1D("phi","#phi","",NBINS,-TMath::Pi(),TMath::Pi());
  // hadron kinematics
  DefineHist1D("z","z","",NBINS,0,1);
  DefineHist1D("pT","p_{T}","GeV",NBINS,1e-4,10);
  DefineHist1D("qT","q_{T}","GeV",NBINS,1e-4,10);
  DefineHist1D("qTq","q_{T}/Q","",NBINS,1e-4,10);
  DefineHist1D("mX","m_{X}","GeV",NBINS,0,20);
  DefineHist1D("phiH","#phi_{h}","",NBINS,-TMath::Pi(),TMath::Pi());
  DefineHist1D("phiS","#phi_{S}","",NBINS,-TMath::Pi(),TMath::Pi());
  // ===========================================================
};


// define a 1D histogram
void Histos::DefineHist1D(
    TString varname,
    TString vartitle,
    TString units,
    Int_t numBins, Double_t lowerBound, Double_t upperBound,
    Bool_t logx
    ) {
  if(units!="") units=" ["+units+"]";
  TH1D *hist = new TH1D(
      setname+"_hist_"+varname,
      vartitle+" distribution, "+settitle+";"+vartitle+units,
      numBins,lowerBound,upperBound
      );
  if(logx) BinLog(hist->GetXaxis());
  histList->AddLast(hist);
  histMap.insert(std::pair<TString,TH1D*>(varname,hist));
  VarNameList.push_back(varname);
};


// define a 2D histogram
void Histos::DefineHist2D(
    TString varname,
    TString vartitlex, TString vartitley,
    TString unitsx, TString unitsy,
    Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
    Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
    Bool_t logx,
    Bool_t logy
    ) {
  if(unitsx!="") unitsx=" ["+unitsx+"]";
  if(unitsy!="") unitsy=" ["+unitsy+"]";
  TH2D *hist = new TH2D(
      setname+"_hist_"+varname,
      " "+vartitley+" vs. "+vartitlex+" distribution, "+settitle+
      ";"+vartitlex+unitsx+";"+vartitley+unitsy,
      numBinsx,lowerBoundx,upperBoundx,
      numBinsy,lowerBoundy,upperBoundy
      );
  if(logx) BinLog(hist->GetXaxis());
  if(logy) BinLog(hist->GetYaxis());
  histList->AddLast(hist);
  histMap.insert(std::pair<TString,TH2D*>(varname,hist));
  VarNameList.push_back(varname);
};


// access histogram by name
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


// make equal-width log-scale bins
void Histos::BinLog(TAxis *axis) {
  Float_t lb = axis->GetXmin();
  Float_t ub = axis->GetXmax();
  if(lb<=0||ub<=0||lb>=ub) {
    fprintf(stderr,"ERROR: bad axis range for Tools::BinLog\n");
    return;
  };
  lb = TMath::Log10(lb);
  ub = TMath::Log10(ub);
  Int_t nb = axis->GetNbins();
  Float_t wd = (ub-lb)/nb;
  Float_t * newBins = new Float_t[nb+1];
  for (int b=0; b<=nb; b++) newBins[b] = TMath::Power(10,lb+b*wd);
  axis->Set(nb,newBins);
  delete[] newBins;
}; 



Histos::~Histos() {
  if(histList) delete histList;
};

