#ifndef Analysis_
#define Analysis_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
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

    // number of bins
    const Int_t NBINS = 50;

    // access bin scheme by name
    BinSet *BinScheme(TString varname);

    // add a new bin scheme
    void AddBinScheme(TString varname, TString vartitle);

    // add a new final state bin
    void AddFinalState(TString finalStateN, TString finalStateT, Int_t pid_);

    // add DIS reconstruction method bin
    void AddRecMethod(TString recMethodN, TString recMethodT);

  
    // get Histos object name and title
    TString GetHistosName(int cpt, int cx, int cz, int cq, int cy, int cfs);
    TString GetHistosTitle(int cpt, int cx, int cz, int cq, int cy, int cfs);
    TString GetHistosNameCoverage(int cp, int cx, int ceta, int cq, int cy, int cfs);
    TString GetHistosTitleCoverage(int cp, int cx, int ceta, int cq, int cy, int cfs);
    TString GetHistosNameJets(int cpt, int cz, int cx, int cq, int cy);
    TString GetHistosTitleJets(int cpt, int cz, int cx, int cq, int cy);
    TString GetHistosNameBreitJets(int cpt, int cz, int cx, int cq, int cy, int crec);
    TString GetHistosTitleBreitJets(int cpt, int cz, int cx, int cq, int cy, int crec);

    // if these are true, only take 'diagonal' elements of the multi
    // dimensional array of possible bins; this is useful
    // if you want to check specific bins
    // TODO: generalize
    Bool_t diagonalPtXZ;
    Bool_t diagonalXZQ;

    // additional settings
    Bool_t writeSimpleTree; // if true, write SimpleTree (not binned)
    Long64_t maxEvents; /* default=0, which runs all events;
                         * if > 0, run a maximum number of `maxEvents` events (useful for quick tests)
                         */

    // perform the analysis
    void Execute();

    // set weights
    void SetWeights(Weights const* w) { weight = w; }
    void SetWeightsJet(Weights const* w) { weightJet = w; }

    // tools
    Bool_t CheckDiagonal(int cpt, int cx, int cz, int cq);

  protected:
    void CheckBins(BinSet *bs, std::vector<int> &v, Double_t var);

  private:
    Histos *HS;
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    TString infileName,outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad
    std::map<TString,BinSet*> binSchemes;

    Weights const* weight;
    Weights const* weightJet;

    std::map<int,int> PIDtoEnum;
    std::map<int,TString> finalStateName;
    std::map<int, TString> recMethodName;
  ClassDef(Analysis,1);
};

#endif
