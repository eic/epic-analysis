// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks, Connor Pecar, Duane Byer, Sanghwa Park, Matthew McEneaney, Brian Page

#ifndef Analysis_
#define Analysis_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <functional>
#include <fmt/format.h>

// root
#include "TChain.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TRegexp.h"

// adage
#include "adage/BinSet.h"

// sidis-eic
#include "Histos.h"
#include "HistosDAG.h"
#include "Kinematics.h"
#include "SimpleTree.h"
#include "Weights.h"
#include "CommonConstants.h"

class Analysis : public TNamed
{
  public:
    Analysis(
        TString configFileName_="",
        TString outfilePrefix_=""
        );
    ~Analysis();

    // number of bins for histograms
    const Int_t NBINS = 50;
    const Int_t NBINS_FULL = 10;

    // bin schemes
    void AddBinScheme(TString varname); // add a new bin scheme
    BinSet *BinScheme(TString varname); // access bin scheme by name
    std::map<TString,BinSet*> GetBinSchemes() { return binSchemes; }; // get full set of bin schemes

    // add a new final state bin
    void AddFinalState(TString finalStateN);

    // common settings
    Bool_t verbose; // if true, print a lot more information
    Bool_t writeSimpleTree; // if true, write SimpleTree (not binned)
    Long64_t maxEvents; /* default=0, which runs all events;
                         * if > 0, run a maximum number of `maxEvents` events (useful for quick tests)
                         */
    Bool_t useBreitJets; // if true, use Breit jets, if using finalState `jets` (requires centauro)
    // set kinematics reconstruction method; see constructor for available methods
    void SetReconMethod(TString reconMethod_) { reconMethod=reconMethod_; }; 
    // choose which output sets to include
    std::map<TString,Bool_t> includeOutputSet;
    // maximum number of errors to print
    Long64_t errorCntMax;

    // Jet Definition Quantities
    int jetAlg;
    double jetRad, jetMin;
    double jetMatchDR;

    // add a group of files to the analysis, where all of these files have a
    // common cross section `xs`, and Q2 range `Q2min` to `Q2max`
    void AddFileGroup(
        std::vector<std::string> fileNames,
        std::vector<Long64_t> entries,
        Double_t xs,
        Double_t Q2min,
        Double_t Q2max
        );

    // access HistosDAG
    HistosDAG *GetHistosDAG() { return HD; };


    // set weights // TODO: are these used yet?
    // void SetWeightsTrack(Weights const* w) { weightTrack = w; }
    // void SetWeightsJet(Weights const* w) { weightJet = w; }

    Double_t GetEventQ2Weight(Double_t Q2, Int_t guess=0);
    // after adding all files, estimate the weights for events in each Q2 range
    void CalculateEventQ2Weights();
    Int_t GetEventQ2Idx(Double_t Q2, Int_t guess=0);

    // run the analysis
    virtual void Execute() = 0;

  protected:

    // prepare to perform the analysis; in derived classes, define a method `Execute()`, which
    // will run the event loop; the first line of `Execute()` should call `Analysis::Prepare()`,
    // which set up common things like output files, `HistosDAG`, etc.
    void Prepare();

    // finish the analysis; call `Analysis::Finish()` at the end of derived `Execute()` methods
    void Finish();

    // print an error; if more than `errorCntMax` errors are printed, printing is suppressed
    void ErrorPrint(std::string message);

    // `FillHistos(weight)` methods: fill histograms
    void FillHistosInclusive(Double_t wgt); // inclusive kinematics
    void FillHistos1h(Double_t wgt);        // single-hadron kinematics
    void FillHistosJets(Double_t wgt);      // jet kinematics

    // shared objects
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    HistosDAG *HD;
    Weights const* weightInclusive;
    Weights const* weightTrack;
    Weights const* weightJet;
    Double_t wInclusiveTotal, wTrackTotal, wJetTotal;
    Long64_t entriesTot;
    Long64_t errorCnt;
    const TString sep = "--------------------------------------------";

    // setup / common settings
    std::vector<std::vector<std::string> > infiles;
    std::vector<std::vector<Long64_t> > inEntries;
    // A lookup index for guessing which Q2 range an event belongs to.
    std::vector<std::size_t> inLookup;
    std::vector<Double_t> Q2xsecs;
    std::vector<Double_t> Q2mins;
    std::vector<Double_t> Q2maxs;
    std::vector<Long64_t> Q2entries;
    std::vector<Double_t> Q2weights;
    TString configFileName,outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn; // GeV
    Double_t ionBeamEn; // GeV
    Double_t crossingAngle; // mrad
    Double_t totalCrossSection;
    TString reconMethod;

    // event loop objects
    Long64_t ENT;
    Double_t eleP,maxEleP;
    Double_t elePtrue, maxElePtrue;
    int pid;
    TString finalStateID;
#ifndef EXCLUDE_DELPHES
    fastjet::PseudoJet jet;
#endif

    // binning names / titles / etc.
    std::map<TString,TString> availableBinSchemes;
    std::map<TString,BinSet*> binSchemes;
    std::map<TString,TString> reconMethodToTitle;
    std::map<TString,TString> finalStateToTitle;
    std::map<int,TString> PIDtoFinalState;
    std::set<TString> activeFinalStates;

    // check if Q2 `val` is between `min` and `max`; if `max==0`, only `val>=min` is checked
    template<class T> bool InQ2Range(T val, T min, T max, bool ignoreZero=false) {
      if (ignoreZero && !(val>0)) return true;
      if (max>0.0) return val>=min && val<=max;
      else         return val>=min;
    }

    // container printing
    // mostly for debugging; if we need more than this, switch to using a common pretty printer library
    template<class O> void PrintStdVector(std::vector<O> vec, std::string name="") {
      if(name!="") fmt::print("{}: ",name);
      fmt::print("[ ");
      for(const auto elem : vec) fmt::print("{}, ",elem);
      fmt::print("]\n");
    }
    template<class O> void PrintStdVector2D(std::vector<std::vector<O>> vec, std::string name="") {
      if(name!="") fmt::print("{} = ",name);
      fmt::print("[\n");
      for(const auto elem : vec) PrintStdVector(elem,"  ");
      fmt::print("]\n");
    }
    template<class K, class V> void PrintStdMap(std::map<K,V> hash, std::string name="") {
      if(name!="") fmt::print("{} = ",name);
      fmt::print("{{\n");
      for(const auto [key,val] : hash) fmt::print("  {} => {},\n",key,val);
      fmt::print("}}\n");
    }

  private:

    // fill histograms, according to `fill_payload`
    void FillHistos(std::function<void(Histos*)> fill_payload);

  ClassDef(Analysis,1);
};

#endif
