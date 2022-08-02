#ifndef BruAsymmetry_
#define BruAsymmetry_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <map>
#include <thread>

// ROOT
#include <TSystem.h>
#include <TStyle.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>

// RooFit
#include <RooDataSet.h>
#include <RooCategory.h>

// BruFit
#include <FitManager.h>
#include <RooMcmc.h>
#include <Process.h>


// DiSpin
#include "Constants.h"
#include "Tools.h"
#include "Modulation.h"


class BruAsymmetry : public TObject
{
  public:
    BruAsymmetry(TString outdir_, TString minimizer_);
    ~BruAsymmetry();

    void AddNumerMod(Modulation * modu);
    void AddDenomMod(Modulation * modu);
    void BuildPDF();
    void LoadDataSets(
        TString dataFileN,
        TString mcFileN="",
        TString weightFileN="",
        TString weightName="Weight",
        TString treeName="tree"
        );
    void Bin(TString varName, Int_t nBins, Double_t *binsArray);
    void Fit();
    void PrintFitter() { FM->SetUp().WS().Print("v"); };

    Int_t GetNdim();
    Int_t GetNbins();
    void PrintBinScheme();

    void PrintLog(TString logString) {
      gSystem->RedirectOutput(outlog,"a");
      printf("%s\n",logString.Data());
      gSystem->RedirectOutput(0);
    };
    TString GetLogName() { return outlog; };

    // hyperparameters
    // - MCMC settings
    Int_t MCMC_iter; // number of MCMC MH steps
    Int_t MCMC_burnin; // discard the first `MCMC_burnin` steps ("burn-in")
    Float_t MCMC_norm; // ~ 1/stepSize
    // - 2nd MCMC settings, with covariance from 1st chain (for MCMCthenCov algorithm)
    Int_t MCMC_cov_iter; // number of MCMC MH steps
    Int_t MCMC_cov_burnin; // discard the first `MCMC_burnin` steps ("burn-in")
    Float_t MCMC_cov_norm; // ~ 1/stepSize

    
    HS::FIT::FitManager * FM;

  private:

    TString outdir;
    TString minimizer;
    TString outlog;

    TFile * infile[2];
    TFile * outfile[2];
    TTree * intr[2];
    TTree * outtr[2];
    Double_t Idx[2];
    TBranch * IdxBr[2];

    TString numerList, PDFstr;
    Int_t nThreads, nWorkers;
    TString numerFormu, denomFormu, formu;
    TString ampNameList,formuNameList;
    Int_t nDenomParams;

    Bool_t useMCint,useWeights;


  ClassDef(BruAsymmetry,1);
};

#endif
