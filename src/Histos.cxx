#include "Histos.h"

ClassImp(Histos)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
Histos::Histos(TString setname_, TString settitle_)
  : setname(setname_)
  , settitle(settitle_)
{
  this->SetName(setname);
  if(settitle!="settitle") cout << "Histos:  " << settitle << endl;
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
  if(logx) BinSet::BinLog(hist->GetXaxis());
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
  if(logx) BinSet::BinLog(hist->GetXaxis());
  if(logy) BinSet::BinLog(hist->GetYaxis());
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
         << histName << " histogram" << endl;
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

// get a specific CutDef
CutDef *Histos::GetCutDef(TString varName) {
  for(auto cut : CutDefList) {
    if(cut->GetVarName() == varName) return cut;
  };
  cerr << "ERROR: cannot find cut " << varName << " in Histos" << endl;
  return nullptr;
};




Histos::~Histos() {
};

