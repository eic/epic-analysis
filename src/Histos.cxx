#include "Histos.h"

ClassImp(Histos)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

Histos::Histos(TString setname_, TString settitle_) {
  setname=setname_;
  settitle=settitle_;
  this->SetName("histos_"+setname);
  if(settitle!="settitle") cout << "Histos:  " << settitle << endl;

  // HISTOGRAMS ================================================
  // -- DIS kinematics
  DefineHist2D("Q2vsX","x","Q^{2}","","GeV^{2}",
      NBINS,1e-3,1,
      NBINS,1,100,
      true,true
      );
  DefineHist1D("W","W","GeV",NBINS,0,15);
  DefineHist1D("y","y","",NBINS,1e-5,1,true);
  //DefineHist1D("Q","Q","GeV",10,0.5,10.5,false,true); // linear
  DefineHist1D("Q","Q","GeV",10,1.0,10.0,true,true); // log
  // -- hadron 4-momentum
  DefineHist1D("p","p","GeV",NBINS,0,10);
  DefineHist1D("pTlab","p_{T}^{lab}","GeV",NBINS,1e-2,3,true);
  DefineHist1D("eta","#eta","",NBINS,-5,5);
  DefineHist1D("phi","#phi","",NBINS,-TMath::Pi(),TMath::Pi());
  // -- hadron kinematics
  DefineHist1D("z","z","",NBINS,0,1);
  DefineHist1D("pT","p_{T}","GeV",NBINS,1e-2,3,true);
  DefineHist1D("qT","q_{T}","GeV",NBINS,1e-2,5,true);
  DefineHist1D("qTq","q_{T}/Q","",NBINS,1e-2,3,true);
  DefineHist1D("mX","m_{X}","GeV",NBINS,0,20);
  DefineHist1D("phiH","#phi_{h}","",NBINS,-TMath::Pi(),TMath::Pi());
  DefineHist1D("phiS","#phi_{S}","",NBINS,-TMath::Pi(),TMath::Pi());
  // -- cross sections
  //DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
  DefineHist1D("Q_xsec","Q","GeV",10,1.0,10.0,true,true); // log
  this->Hist("Q_xsec")->SetMinimum(1e-10);
  // ===========================================================
};


// define a 1D histogram
void Histos::DefineHist1D(
    TString varname,
    TString vartitle,
    TString units,
    Int_t numBins, Double_t lowerBound, Double_t upperBound,
    Bool_t logx, Bool_t logy
    ) {
  if(units!="") units=" ["+units+"]";
  TString histT;
  if(varname.Contains("_xsec")) histT = "d#sigma/d"+vartitle;
  else if(varname.Contains("_fuu")) histT = "F_{UU} vs. "+vartitle;
  else if(varname.Contains("_fut")) histT = "F_{UT} vs. "+vartitle;
  else histT = vartitle+" distribution";
  TH1D *hist = new TH1D(
      setname+"_hist_"+varname,
      histT+", "+settitle+";"+vartitle+units,
      numBins,lowerBound,upperBound
      );
  if(logx) BinLog(hist->GetXaxis());
  HistConfig *config = new HistConfig();
  config->logx = logx;
  config->logy = logy;
  config->logz = false;
  this->RegisterHist(varname,hist,config);
};


// define a 2D histogram
void Histos::DefineHist2D(
    TString varname,
    TString vartitlex, TString vartitley,
    TString unitsx, TString unitsy,
    Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
    Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
    Bool_t logx, Bool_t logy, Bool_t logz
    ) {
  if(unitsx!="") unitsx=" ["+unitsx+"]";
  if(unitsy!="") unitsy=" ["+unitsy+"]";
  TString histT;
  if(varname.Contains("_xsec")) histT = "d^{2}#sigma/d"+vartitlex+vartitley;
  else if(varname.Contains("_fuu")) histT = "F_{UU} vs. "+vartitley+" vs. "+vartitlex;
  else if(varname.Contains("_fut")) histT = "F_{UT} vs. "+vartitley+" vs. "+vartitlex;
  else histT = vartitley+" vs. "+vartitlex+" distribution";
  TH2D *hist = new TH2D(
      setname+"_hist_"+varname,
      histT+", "+settitle+";"+vartitlex+unitsx+";"+vartitley+unitsy,
      numBinsx,lowerBoundx,upperBoundx,
      numBinsy,lowerBoundy,upperBoundy
      );
  if(logx) BinLog(hist->GetXaxis());
  if(logy) BinLog(hist->GetYaxis());
  HistConfig *config = new HistConfig();
  config->logx = logx;
  config->logy = logy;
  config->logz = logz;
  this->RegisterHist(varname,hist,config);
};


// add histogram to containers
void Histos::RegisterHist(TString varname_, TH1 *hist_, HistConfig *config_) {
  VarNameList.push_back(varname_);
  histConfigMap.insert(std::pair<TString,HistConfig*>(varname_,config_));
  histMap.insert(std::pair<TString,TH1*>(varname_,hist_));
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
// access histogram config by name
HistConfig *Histos::GetHistConfig(TString histName) {
  HistConfig *retConfig;
  try { retConfig = histConfigMap.at(histName); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: histConfigMap does not have " 
         << histName << "histogram" << endl;
    return nullptr;
  };
  return retConfig;
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
};

