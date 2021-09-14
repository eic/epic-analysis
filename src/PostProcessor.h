#ifndef PostProcessor_
#define PostProcessor_

#include <map>
#include <iomanip>

// root
#include "TFile.h"
#include "TCanvas.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"

// largex-eic
#include "Histos.h"
#include "Kinematics.h"
#include "CutDef.h"
#include "BinSet.h"
#include "HistosDAG.h"

class PostProcessor : public TNamed
{
  public:
    PostProcessor(
        TString infileN_
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
    void DrawSingle(Histos *H, TString histName, TString drawFormat="");
    void DrawSingle(TString histSet, TString histName);
    void DrawRatios(
        TString outName, Histos *numerSet, Histos *denomSet, Bool_t plotRatioOnly=false
        );
    void DrawInBins(
        TString outName,
        std::vector<std::vector<Histos*>>& histList, TString histName,
        TString var1name, int nvar1, double var1low, double var1high, bool var1log,
        TString var2name, int nvar2, double var2low, double var2high, bool var2log
        );

    // algorithm finish methods; to be called after loops
    void FinishDumpAve(TString datFile);
    void FinishDrawRatios(TString summaryDir);

    // accessors
    TString GetPngDir() { return pngDir; };
    TString GetOutfileName() { return outfileN; };
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

#endif
