#ifndef Analysis_
#define Analysis_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <stdexcept>

// root
#include "TChain.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TRegexp.h"

// delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

//#include "fastjet/contrib/Centauro.hh"
//#include "fastjet/plugins/Centauro/Centauro.hh"

// largex-eic
#include "Histos.h"
#include "HistosDAG.h"
#include "Kinematics.h"
#include "CutDef.h"
#include "BinSet.h"
#include "SimpleTree.h"
#include "Weights.h"


class Analysis : public TNamed
{
  public:
    Analysis(
        TString infileName_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0,
        TString outfilePrefix_=""
        );
    ~Analysis();

    // number of bins for histograms
    const Int_t NBINS = 50;

    // bin schemes
    void AddBinScheme(TString varname); // add a new bin scheme
    BinSet *BinScheme(TString varname); // access bin scheme by name
    std::map<TString,BinSet*> GetBinSchemes() { return binSchemes; }; // get full set of bin schemes

    // add a new final state bin
    void AddFinalState(TString finalStateN);

    // additional settings
    Bool_t writeSimpleTree; // if true, write SimpleTree (not binned)
    Long64_t maxEvents; /* default=0, which runs all events;
                         * if > 0, run a maximum number of `maxEvents` events (useful for quick tests)
                         */
    Bool_t useBreitJets; // if true, use Breit jets, if using finalState `jets` (requires centauro)
    // set kinematics reconstruction method; see constructor for available methods
    void SetReconMethod(TString reconMethod_) { reconMethod=reconMethod_; }; 

    // perform the analysis
    void Execute();

    // set weights
    void SetWeights(Weights const* w) { weight = w; }
    void SetWeightsJet(Weights const* w) { weightJet = w; }

  protected:
    void CheckBins(BinSet *bs, std::vector<int> &v, Double_t var);

  private:

    Histos *HS;
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    HistosDAG *HD;

    TString infileName,outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad

    TString reconMethod;
    Double_t eleP,maxEleP;
    Double_t elePtrue, maxElePtrue;
    int pid;

    std::map<TString,BinSet*> binSchemes;
    std::map<TString,TString> availableBinSchemes;
    std::map<TString,TString> reconMethodToTitle;

    std::map<TString,Double_t> valueMap;
    TString finalStateID;
    Bool_t activeEvent;

    Weights const* weight;
    Weights const* weightJet;

    std::map<TString, TString> finalStateToTitle;
    std::map<int, TString> PIDtoFinalState;
    std::set<TString> activeFinalStates;
  ClassDef(Analysis,1);
};

#endif
