#ifndef Analysis_
#define Analysis_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <stdexcept>
#include <functional>

// root
#include "TChain.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TRegexp.h"

// largex-eic
#include "Histos.h"
#include "HistosDAG.h"
#include "Kinematics.h"
#include "CutDef.h"
#include "BinSet.h"
#include "SimpleTree.h"
#include "Weights.h"

// delphes (TODO: does fastjet need this?)
//#include "classes/DelphesClasses.h"
//#include "external/ExRootAnalysis/ExRootTreeReader.h"



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

    // common settings
    Bool_t writeSimpleTree; // if true, write SimpleTree (not binned)
    Long64_t maxEvents; /* default=0, which runs all events;
                         * if > 0, run a maximum number of `maxEvents` events (useful for quick tests)
                         */
    Bool_t useBreitJets; // if true, use Breit jets, if using finalState `jets` (requires centauro)
    // set kinematics reconstruction method; see constructor for available methods
    void SetReconMethod(TString reconMethod_) { reconMethod=reconMethod_; }; 

    // add files to the TChain; this is called by `Prepare()`, but you can use these public
    // methods to add more files if you want
    void AddFile(TString fileName); // add single file `fileName`
    void AddFiles(TString fileList); // add files listed in `fileList`

    // access HistosDAG
    HistosDAG *GetHistosDAG();

    // set weights // TODO: are these used yet?
    void SetWeights(Weights const* w) { weight = w; }
    void SetWeightsJet(Weights const* w) { weightJet = w; }

    // run the analysis
    virtual void Execute() = 0;

  protected:

    // prepare to perform the analysis; in derived classes, define a method `Execute()`, which
    // will run the event loop; the first line of `Execute()` should call `Analysis::Prepare()`,
    // which set up common things like output files, `HistosDAG`, etc.
    void Prepare();

    // finish the analysis; call `Analysis::Finish()` at the end of derived `Execute()` methods
    void Finish();

    // FillHistos methods: fill histograms
    void FillHistosTracks();
    void FillHistosJets();

    // lambda to check which bins an observable is in, during DAG breadth
    // traversal; it requires `finalStateID`, `valueMap`, and will
    // activate/deactivate bin nodes accoding to values in `valuMap`
    std::function<void(Node*)> CheckBin();


    // shared object pointers
    Histos *HS;
    SimpleTree *ST;
    Kinematics *kin, *kinTrue;
    HistosDAG *HD;
    Weights const* weight;
    Weights const* weightJet;

    // setup / common settings
    std::vector<TString> infiles;
    TString infileName,outfileName,outfilePrefix;
    TFile *outFile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad
    TString reconMethod;

    // event loop objects
    Long64_t ENT;
    Double_t eleP,maxEleP;
    Double_t elePtrue, maxElePtrue;
    int pid;
    fastjet::PseudoJet jet;
    std::map<TString,Double_t> valueMap;
    TString finalStateID;
    Bool_t activeEvent;
    Double_t wTrack,wJet;

    // binning names / titles / etc.
    std::map<TString,TString> availableBinSchemes;
    std::map<TString,BinSet*> binSchemes;
    std::map<TString,TString> reconMethodToTitle;
    std::map<TString, TString> finalStateToTitle;
    std::map<int, TString> PIDtoFinalState;
    std::set<TString> activeFinalStates;

  ClassDef(Analysis,1);
};

#endif
