// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Duane Byer

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
// define 2D histogram with custom bins
void Histos::DefineHist2D(
    TString varname,
    TString vartitlex, TString vartitley,
    TString unitsx, TString unitsy,
    Int_t numBinsx, Double_t *xBins, 
    Int_t numBinsy, Double_t *yBins, 
    Bool_t logx, Bool_t logy, Bool_t logz
    ) {
  if(unitsx!="") unitsx=" ["+unitsx+"]";
  if(unitsy!="") unitsy=" ["+unitsy+"]";
  TString histT;
  TH2D *hist = new TH2D(
      setname+"_hist_"+varname,
      histT+", "+settitle+";"+vartitlex+unitsx+";"+vartitley+unitsy,
      numBinsx,xBins,
      numBinsy,yBins
      );
  HistConfig *config = new HistConfig();
  config->logx = logx;
  config->logy = logy;
  config->logz = logz;
  this->RegisterHist(varname,hist,config);
};

// define a 3D histogram
void Histos::DefineHist3D(
    TString varname,
    TString vartitlex, TString vartitley, TString vartitlez,
    TString unitsx, TString unitsy, TString unitsz,
    Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
    Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
    Int_t numBinsz, Double_t lowerBoundz, Double_t upperBoundz,
    Bool_t logx, Bool_t logy, Bool_t logz
    ) {
  if(unitsx!="") unitsx=" ["+unitsx+"]";
  if(unitsy!="") unitsy=" ["+unitsy+"]";
  if(unitsz!="") unitsz=" ["+unitsz+"]";
  TString histT;
  if(varname.Contains("_xsec")) histT = "d^{3}#sigma/d"+vartitlex+vartitley+vartitlez;
  else if(varname.Contains("_fuu")) histT = "F_{UU} vs. "+vartitlez+" vs. "+vartitley+" vs. "+vartitlex;
  else if(varname.Contains("_fut")) histT = "F_{UT} vs. "+vartitlez+" vs. "+vartitley+" vs. "+vartitlez;
  else histT = vartitlez+" vs. "+vartitley+" vs. "+vartitlex+" distribution";
  TH3D *hist = new TH3D(
      setname+"_hist_"+varname,
      histT+", "+settitle+";"+vartitlex+unitsx+";"+vartitley+unitsy+";"+vartitlez+unitsz,
      numBinsx,lowerBoundx,upperBoundx,
      numBinsy,lowerBoundy,upperBoundy,
      numBinsz,lowerBoundz,upperBoundz
      );
  if(logx) BinSet::BinLog(hist->GetXaxis());
  if(logy) BinSet::BinLog(hist->GetYaxis());
  if(logz) BinSet::BinLog(hist->GetZaxis());
  HistConfig *config = new HistConfig();
  config->logx = logx;
  config->logy = logy;
  config->logz = logz;
  this->RegisterHist(varname,hist,config);
};


// define a 4D histogram
void Histos::DefineHist4D(
    TString varname,
    TString vartitlew, TString vartitlex, TString vartitley, TString vartitlez,
    TString unitsw, TString unitsx, TString unitsy, TString unitsz,
    Int_t numBinsw, Double_t lowerBoundw, Double_t upperBoundw,
    Int_t numBinsx, Double_t lowerBoundx, Double_t upperBoundx,
    Int_t numBinsy, Double_t lowerBoundy, Double_t upperBoundy,
    Int_t numBinsz, Double_t lowerBoundz, Double_t upperBoundz,
    Bool_t logw, Bool_t logx, Bool_t logy, Bool_t logz
    ) {
  if(unitsw!="") unitsw=" ["+unitsw+"]";
  if(unitsx!="") unitsx=" ["+unitsx+"]";
  if(unitsy!="") unitsy=" ["+unitsy+"]";
  if(unitsz!="") unitsz=" ["+unitsz+"]";
  TString histT;
  if(varname.Contains("_xsec")) histT = "d^{4}#sigma/d"+vartitlew+vartitlex+vartitley+vartitlez;
  else if(varname.Contains("_fuu")) histT = "F_{UU} vs. "+vartitlez+" vs. "+vartitley+" vs. "+vartitlex+" vs. "+vartitlew;
  else if(varname.Contains("_fut")) histT = "F_{UT} vs. "+vartitlez+" vs. "+vartitley+" vs. "+vartitlez+" vs. "+vartitlew;
  else histT = vartitlez+" vs. "+vartitley+" vs. "+vartitlex+" vs. "+vartitlew+" distribution";
  Hist4D *hist = new Hist4D(
      setname+"_hist_"+varname,
      histT+", "+settitle+";"+vartitlew+unitsw+";"+vartitlex+unitsx+";"+vartitley+unitsy+";"+vartitlez+unitsz,
      numBinsw,lowerBoundw,upperBoundw,
      numBinsx,lowerBoundx,upperBoundx,
      numBinsy,lowerBoundy,upperBoundy,
      numBinsz,lowerBoundz,upperBoundz
      );
  if(logw) BinSet::BinLog(hist->GetWaxis());
  if(logx) BinSet::BinLog(hist->GetXaxis());
  if(logy) BinSet::BinLog(hist->GetYaxis());
  if(logz) BinSet::BinLog(hist->GetZaxis());
  HistConfig *config = new HistConfig();
  config->logw = logw;
  config->logx = logx;
  config->logy = logy;
  config->logz = logz;
  this->RegisterHist4(varname,hist,config);
};


// add histogram to containers
void Histos::RegisterHist(TString varname_, TH1 *hist_, HistConfig *config_) {
  VarNameList.push_back(varname_);
  histConfigMap.insert({varname_,config_});
  histMap.insert({varname_,hist_});
};

void Histos::RegisterHist4(TString varname_, Hist4D *hist_, HistConfig *config_) {
  VarNameList.push_back(varname_);
  hist4ConfigMap.insert({varname_,config_});
  hist4Map.insert({varname_,hist_});
};


// access histogram by name
TH1 *Histos::Hist(TString histName, Bool_t silence) {
  TH1 *retHist;
  try { retHist = histMap.at(histName); }
  catch(const std::out_of_range &ex) {
    if(!silence)
      cerr << "ERROR: histMap does not have " 
           << histName << " histogram" << endl;
    return nullptr;
  };
  return retHist;
};

Hist4D *Histos::Hist4(TString histName, Bool_t silence) {
  Hist4D *retHist;
  try { retHist = hist4Map.at(histName); }
  catch(const std::out_of_range &ex) {
    if(!silence)
      cerr << "ERROR: hist4Map does not have " 
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

HistConfig *Histos::GetHist4Config(TString histName) {
  HistConfig *retConfig;
  try { retConfig = hist4ConfigMap.at(histName); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: hist4ConfigMap does not have " 
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

// add histograms of `in_histos` to histograms in `this`
void Histos::AddHistos(Histos *in_histos) {
  for(auto [key,this_hist] : histMap) {
    auto in_hist = in_histos->Hist(key);
    if(in_hist==nullptr) continue;
    this_hist->Add(in_hist);
  }
  // for(auto [key,this_hist] : hist4Map) {
  //   auto in_hist = in_histos->Hist4(key);
  //   if(in_hist==nullptr) continue;
  //   this_hist->Add(in_hist); // TODO: no add method
  // }
}

Histos::~Histos() {
  for(auto it : CutDefList)
    if(it) delete it;
  auto del = [] (auto m) {
    for(auto [k,v] : m)
      if(v) delete v;
    m.clear();
  };
  del(histMap);
  del(hist4Map);
  del(histConfigMap);
  del(hist4ConfigMap);
}
