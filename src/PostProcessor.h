// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Connor Pecar, Duane Byer

#pragma once

#include <map>
#include <iomanip>

// root
#include <TFile.h>
#include <TCanvas.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TProfile2D.h>

// adage
#include <adage/CutDef.h>
#include <adage/BinSet.h>

// epic-analysis
#include "Histos.h"
#include "Kinematics.h"
#include "HistosDAG.h"

class PostProcessor
{
  public:
    PostProcessor(
        TString infileN_,
        TString outfileN_=""
        );
    ~PostProcessor();

    // settings
    // - canvas dimensions [pixels]
    const Int_t dimx=800;
    const Int_t dimy=700;
    static const int nsumMax=3; // number of summary plots with formatting


    // DAG interfaces:
    HistosDAG *GetHistosDAG() { return HD; };
    HistosDAG *Op() { return GetHistosDAG(); }; // syntactic sugar
    // execute lambdas (if `clear`==false, lambda operators will not be removed after execution)
    void Execute(Bool_t clear=true) {
      if(clear) HD->ExecuteAndClearOps();
      else HD->ExecuteOps();
    };


    // cleanup and close open files and streams
    // - MUST be called at the end of any postprocessor macro, after
    //   all algorithms have finished
    void Finish();


    // algorithms: useful to run in loops over bins
    // - see PostProcessor.cxx for descriptions for how to use these
    // - these are general functions that operate either on Histos 
    //   objects, or on histograms
    // - they can be shared in any postprocessor macro
    // - they can be anything, such as taking ratios of each histogram
    //   from two different Histos objects (e.g., y>0.05 / y>0.00 sets)
    // - you are welcome to add your own algorithms
    void DumpHist(TString datFile, TString histSet, TString varName);
    void DumpAve(TString datFile, Histos *H, TString cutName);
    void DrawSingle(Histos *H, TString histName, TString drawFormat="", Int_t profileAxis=0, Bool_t profileOnly=false);
    void DrawSingle(TString histSet, TString histName);
    void DrawRatios(
        TString outName, Histos *numerSet, Histos *denomSet, Bool_t plotRatioOnly=false
        );
    void DrawInBins(
        TString outName,
        std::vector<std::vector<Histos*>>& histArr, TString histName,
        TString var1name, int nvar1, double var1low, double var1high, bool var1log,
        TString var2name, int nvar2, double var2low, double var2high, bool var2log,
        bool intgrid1=false, bool intgrid2=false,
        bool renormalize=false
        );
    void DrawInBins(
        TString outName,
        std::vector<std::vector<std::vector<Histos*>>>& histArrList, TString histName,
        TString var1name, int nvar1, double var1low, double var1high, bool var1log,
        TString var2name, int nvar2, double var2low, double var2high, bool var2log,
        bool intgrid1=false, bool intgrid2=false,
        bool renormalize=false
        );
    void DrawRatioInBins(
        TString outName,
        std::vector<std::vector<Histos*>>& histArr, TString histName,
        TString var1name, int nvar1, double var1low, double var1high, bool var1log,
        TString var2name, int nvar2, double var2low, double var2high, bool var2log,
        bool intgrid1=false, bool intgrid2=false,
        bool renormalize=false
        );
    void DrawRatioInBins(
        TString outName,
        std::vector<std::vector<std::vector<Histos*>>>& histArrList, TString histName,
        TString var1name, int nvar1, double var1low, double var1high, bool var1log,
        TString var2name, int nvar2, double var2low, double var2high, bool var2log,
        bool intgrid1=false, bool intgrid2=false,
        bool renormalize=false
        );

    // algorithm finish methods; to be called after loops
    void FinishDumpAve(TString datFile);
    void FinishDrawRatios(TString summaryDir);

    // vector of labels for a legend (push elements externally if you want to use this)
    std::vector<TString> legendLabels;

    // accessors
    TString GetPngDir() { return pngDir; };
    TString GetOutfileName() { return outfileN; };
    TFile *GetOutfile() { return outfile; };
    BinSet *GetBinSet(TString varName);
    CutDef *GetBinCut(TString varName, Int_t binNum);
    std::vector<int> GetBinNums(TString varName);

    // text file manipulation
    void StartTextFile(TString datFile, TString firstLine="");
    void AppendToTextFile(TString datFile, TString appendText);
    void Columnify(TString inputFile, TString outputFile);
    void PrintTextFile(TString datFile);

    // return true if the bin is "full" range, and it's not the only bin
    Bool_t SkipFull(TString varName, Int_t binNum);

    // reset algorithm-specific variables
    void ResetVars();

    
    // ------------
    // zoom out the vertical scale for the case where multiple
    // `TH1`s have been drawn with the "SAME" option, but the y-axis
    // range is improperly zoomed
    // - example: `UnzoomVertical(canvas->GetPad(3))`
    // - optionally specify a new title 
    // - set `min0` to true if you want to lock the minimum at zero
    static void UnzoomVertical(TVirtualPad *pad, TString title="", Bool_t min0=false, Bool_t logy=false) {
      Double_t max=-1e6;
      Double_t min=1e6;
      Double_t maxTmp,minTmp;
      for(auto obj : *pad->GetListOfPrimitives()) {
        if(obj->InheritsFrom(TH1::Class())) {
          maxTmp = ((TH1*)obj)->GetMaximum();
          minTmp = ((TH1*)obj)->GetMinimum();
          max = maxTmp > max ? maxTmp : max;
          min = minTmp < min ? minTmp : min;
        };
      };
      max += 0.05*(max-min);
      //min -= 0.05*(max-min);
      for(auto obj : *pad->GetListOfPrimitives()) {
        if(obj->InheritsFrom(TH1::Class())) {
          Double_t drawMin;
          if(logy) drawMin = 0.001*max;
          else drawMin = min0 ? 0:min;
          ((TH1*)obj)->GetYaxis()->SetRangeUser(drawMin,max);
          if(title!="") ((TH1*)obj)->SetTitle(title);
        };
      };
    };
    // ------------



  private:

    // files and names
    TString infileN, outfileN, pngDir;
    TFile *infile, *outfile;

    // DAGs
    HistosDAG *HD;

    // algorithm-specific variables
    std::map<TString,TCanvas*> summaryCanvMap;
    std::vector<TString> varList;
    int nsum,ndump;
    CutDef *dumpCut;
    TCanvas *summaryCanv;
    Color_t summaryColor[nsumMax];
    Style_t summaryStyle[nsumMax];

  ClassDef(PostProcessor,1);
};
