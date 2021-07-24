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

// largex-eic
#include "Histos.h"
#include "Kinematics.h"
#include "CutDef.h"
#include "BinSet.h"


class Analysis : public TNamed
{
  public:
    Analysis(
        TString infile_="",
        Double_t eleBeamEn_=5,
        Double_t ionBeamEn_=41,
        Double_t crossingAngle_=0
        );
    ~Analysis();

    // access bin scheme by name
    BinSet *BinScheme(TString varname);

    // add a new bin scheme
    void AddBinScheme(TString varname, TString vartitle);

    // if this is true, only take 'diagonal' elements of the multi
    // dimensional array of possible pT,x,z bins; this is useful
    // if you want to check specific bins
    Bool_t diagonalBinsOnly;

    // perform the analysis
    void Execute();


  protected:
    Bool_t CheckDiagonal(int cpt, int cx, int cz);
    void CheckBins(BinSet *bs, std::vector<int> &v, Double_t var);

  private:
    TString infile;
    Double_t eleBeamEn = 5; // GeV
    Double_t ionBeamEn = 41; // GeV
    Double_t crossingAngle = 0; // mrad
    std::map<TString,BinSet*> binSchemes;
    std::map<int,int> PIDtoEnum;
    enum partEnum{
      pPip,
      //pPim,
      NPart
    };

  ClassDef(Analysis,1);
};

#endif
