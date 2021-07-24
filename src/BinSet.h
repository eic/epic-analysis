#ifndef BinSet_
#define BinSet_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

// ROOT
#include "TSystem.h"
#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TMath.h"
#include "TAxis.h"
#include "TObjArray.h"


// largex-eic
#include "CutDef.h"


class BinSet : public TObject
{
  public:
    /* this class stores `CutDef`s, for your specified binning
     * - `varName` and `varTitle` are passed to the `CutDef` objects
     */
    BinSet(TString varName_="unknown", TString varTitle_="unknown");
    ~BinSet();

    // bin list container, of `CutDef` pointers, one for each bin
    TObjArray *GetBinList() { return binList; };

    /* bin builders
     * - at construction, you will start with zero bins
     * - call any bin builder to sequentiall add bins to the list of bins
     */
    /* build a single bin
     * - bin is created by specifying a `CutDef` (see CutDef.h)
     */
    void BuildBin(TString cutType_, Double_t arg1_=-1, Double_t arg2_=-1);
    void BuildBin(CutDef *cut_);
    /* build list of bins
     * - define the number of bins `nbins`, in the * range `min` to `max`
     * - default is equal width bins in linear scale
     * - set `log` to true for equal-width bins in log scale
     * - you may instead specify a `TAxis`, for any arbitrary binning
     * - a list of `CutDef`s is generated
     */
    void BuildBins(Int_t nbins_, Double_t min_, Double_t max_, Bool_t log_=false);
    void BuildBins(TAxis *ax, Bool_t log_=false);

    // access a bin's cut (by index)
    CutDef *Cut(Int_t idx) { return (CutDef*)binList->At(idx); };

    /* make equal-width log-scale bins
     * - the axis `ax` will be modified
     * - you can use this on histogram axes too; just call this method
     *   on any histogram axis, prior to filling the histogram
     */
    static void BinLog(TAxis *ax);

  private:
    TObjArray *binList;
    TString varName,varTitle;
    Int_t nbins;
    Double_t min,max;
    Bool_t log;

  ClassDef(BinSet,1);
};

#endif
