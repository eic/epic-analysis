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
    /* builds container `bins` of `CutDef`s, for your specified binning with
     * `nbins`, `min`, and `max`
     * - `varName` and `varTitle` are passed to the `CutDef` objects
     * - set `log` to true for equal-width binning in log scale
     */
    BinSet(
        TString varName_, TString varTitle_,
        Int_t nbins_, Double_t min_, Double_t max_, Bool_t log_=false
        );
    ~BinSet();

    /* bin list container, of `CutDef` pointers, one for each bin
     */
    TObjArray *bins;

    /* make equal-width log-scale bins
     * - the axis `ax` will be modified
     * - it is possible to use this on histogram axes; just call this method
     *   on the histogram axis, prior to filling the histogram
     */
    static void BinLog(TAxis *ax);

  private:
    TAxis *axis;
    TString varName,varTitle;
    Int_t nbins;
    Double_t min,max;
    Bool_t log;

  ClassDef(BinSet,1);
};

#endif
