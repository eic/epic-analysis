R__LOAD_LIBRARY(LargexAsym)
#include "BruAsymmetry.h"
// draw brufit results

//////////////////////////////////////////////////////////

TObjArray * BruBinList;
TObjArray * BruBinSuperList;
Int_t nDim, nBins;
Int_t nSamples;
TString bruDir;
HS::FIT::Bins * HSbins;
const int nParamsMax = 30;
Float_t asymPlotMin;
Float_t asymPlotMax;
Float_t axisTitleSize;
Int_t minimizer;
enum minimEnum { mkMCMC, mkMCMCthenCov, mkMinuit };
TString hTitle,vTitle,pName;
TString resultFileN;


//////////////////////////////////////////////////////////

class BruBin : public TObject {
  public:
    Int_t idx;
    TString name,var;
    Double_t center,mean,ub,lb;
    TH1D * hist;
    Double_t param[nParamsMax];
    Double_t paramErr[nParamsMax];
    TH1D * paramVsSample[nParamsMax];
    TH1D * nllVsSample;
    TFile * resultFile;
    TTree * resultTree;
    TTree * mcmcTree;

    // -constructor
    BruBin(TAxis axis, Int_t binnum, TVectorD coord) {

      idx = HSbins->FindBin(coord); // bin index
      name = HSbins->GetBinName(idx); // bin name (brufit syntax)
      var = axis.GetName(); // bin variable

      // get bin center, lower bound, upper bound
      center = axis.GetBinCenter(binnum);
      lb = axis.GetBinLowEdge(binnum);
      ub = axis.GetBinUpEdge(binnum);

      // get bin mean (requires opening the bin's tree)
      hist = new TH1D( Form("ivHist%d",idx),Form("ivHist%d",idx),
        100,lb-0.05,ub+0.05);
      TFile * binTreeFile;
      TTree * binTree;
      Double_t iv;
      binTreeFile = new TFile(bruDir+"/"+name+"/TreeData.root","READ");
      binTree = (TTree*) binTreeFile->Get("tree");
      binTree->SetBranchAddress(var,&iv);
      for(Long64_t e=0; e<binTree->GetEntries(); e++) {
        binTree->GetEntry(e);
        hist->Fill(iv);
      };
      mean = hist->GetMean();
      binTreeFile->Close();

      // open result file, and read trees
      switch(minimizer) {
        case mkMCMC:        resultFileN="ResultsHSRooMcmcSeq.root"; break;
        case mkMCMCthenCov: resultFileN="ResultsHSRooMcmcSeqThenCov.root"; break;
        case mkMinuit:      resultFileN="ResultsHSMinuit2.root"; break;
      };
      resultFile = new TFile(bruDir+"/"+name+"/"+resultFileN,"READ");
      resultTree = (TTree*) resultFile->Get("ResultTree");
      if(minimizer==mkMCMC || minimizer==mkMCMCthenCov) {
        mcmcTree = (TTree*) resultFile->Get("MCMCTree");
        nSamples = (Int_t) mcmcTree->GetEntries();
        for(int i=0; i<nParamsMax; i++) {
          pName = Form("bin%d_param%d_vs_sample",idx,i);
          paramVsSample[i] = new TH1D(pName,pName,nSamples,1,nSamples+1);
        };
        pName = Form("bin%d_NLL_vs_sample",idx);
        nllVsSample = new TH1D(pName,pName,nSamples,1,nSamples+1);
      };

      // misc
      for(int i=0; i<nParamsMax; i++) {
        param[i] = UNDEF;
        paramErr[i] = UNDEF;
      };
    };

    // -print out
    void PrintInfo() {
      printf("BIN %d\n",idx);
      printf("  name   = %s\n",name.Data());
      printf("  var    = %s\n",var.Data());
      printf("  center = %.2f\n",center);
      printf("  mean   = %.2f\n",mean);
      printf("  range  = %.2f to %.2f\n",lb,ub);
    };
};

//////////////////////////////////////////////////////////

void draw(
  TString bruDir_ = "bruspin",
  TString minimizer_ = "minuit",
  Bool_t logscale = true, /* plot sample# on log scale */
  Float_t asymPlotMin_ = -10, /* units: % */
  Float_t asymPlotMax_ =  10  /* units: % */
) {

  asymPlotMin = asymPlotMin_ / 100.0;
  asymPlotMax = asymPlotMax_ / 100.0;
  gStyle->SetTitleSize(0.08,"T");
  axisTitleSize = 0.05;
  gStyle->SetOptStat(0);

  // get minimizer type
  if(minimizer_.CompareTo("mcmc",TString::kIgnoreCase)==0) minimizer=mkMCMC;
  else if(minimizer_.CompareTo("mcmcthencov",TString::kIgnoreCase)==0) minimizer=mkMCMCthenCov;
  else if(minimizer_.CompareTo("minuit",TString::kIgnoreCase)==0) minimizer=mkMinuit;
  else { fprintf(stderr,"ERROR: unknown minimizer type\n"); return; };


  // get binning scheme
  bruDir = bruDir_;
  TFile * binFile = new TFile(bruDir+"/DataBinsConfig.root","READ");
  HSbins = (HS::FIT::Bins*) binFile->Get("HSBins");
  nDim = HSbins->GetNAxis();
  printf("nDim = %d\n",nDim);
  binFile->Close();

  // get Nbins
  nBins = HSbins->GetN();
  printf("nBins = %d\n",nBins);

  // initialize nSamples
  nSamples = 0;

  // build arrays of BruBin objects
  // - BruBinList: bins for the first dimension (IV0)
  // - BruBinSuperList: array of BruBinLists, one for each higher-dimensional bin
  // -- get axes:
  TAxis ax[3];
  int bn[3];
  int bnMax[3];
  for(int i=0; i<3; i++) {
    if(i<nDim) {
      ax[i] = HSbins->GetVarAxis()[i];
      bnMax[i] = ax[i].GetNbins();
    } else bnMax[i] = 1;
  };
  hTitle = ax[0].GetName(); // TODO [low priority]: use dispin Binning to convert to proper title
  // -- create BruBin objects, and fill TObjArrays
  TVectorD binCoord(nDim);
  BruBinSuperList = new TObjArray();
  for(bn[2]=1; bn[2]<=bnMax[2]; bn[2]++) {
    for(bn[1]=1; bn[1]<=bnMax[1]; bn[1]++) {
      BruBinList = new TObjArray(); // create new BruBinList before looping over dimension 0
      for(bn[0]=1; bn[0]<=bnMax[0]; bn[0]++) {
        // get bin center coordinate
        for(int i=0; i<nDim; i++) binCoord[i] = ax[i].GetBinCenter(bn[i]);
        // create new BruBin object
        BruBinList->AddLast(new BruBin(ax[0],bn[0],binCoord));
      }
      BruBinSuperList->AddLast(BruBinList);
    };
  };


  // define BruBin iterator, and print bins
  BruBin * BB;
  TObjArray * BBlist;
  TObjArrayIter nextBin(BruBinList);
  TObjArrayIter nextBinList(BruBinSuperList);
  /*
  Tools::PrintSeparator(30,"-");
  Int_t BLn=0;
  while((BBlist = (TObjArray*) nextBinList())) {
    printf("BINLIST %d\n",BLn++);
    nextBin = TObjArrayIter(BBlist);
    while((BB = (BruBin*) nextBin())) BB->PrintInfo();
    Tools::PrintSeparator(30,"-");
  };
  */


  // get parameter values from results trees
  Bool_t first = true;
  RooDataSet * paramSet;
  TString paramList[nParamsMax];
  Modulation * moduList[nParamsMax];
  Int_t nParams;
  TString paramName;
  Double_t paramval[nParamsMax];
  Double_t nll;
  Long64_t entry;
  nextBinList = TObjArrayIter(BruBinSuperList);
  while((BBlist = (TObjArray*) nextBinList())) {
    nextBin = TObjArrayIter(BBlist);
    while((BB = (BruBin*) nextBin())) {

      // get parameter tree
      paramSet = (RooDataSet*) BB->resultFile->Get("FinalParameters");

      // get parameter list
      if(first) {
        nParams = 0;
        for(int i=0; i<paramSet->get()->size(); i++) {
          if(nParams>nParamsMax) {
            fprintf(stderr,"ERROR: too many params\n"); return 1; };
          paramName = (*(paramSet->get()))[i]->GetName();
          if(paramName=="NLL") continue;
          if(paramName.Contains("Yld")) moduList[nParams] = nullptr;
          else moduList[nParams] = new Modulation(paramName);
          paramList[nParams] = paramName;
          printf("param %d:  %s\n",nParams,paramList[nParams].Data());
          nParams++;
        };
        Tools::PrintSeparator(30);
        first = false;
      };

      // get parameter values
      if(BB->resultTree->GetEntries()!=1)
        fprintf(stderr,"WARNING: ResultTree does not have 1 entry\n");
      for(int i=0; i<nParams; i++) {
        BB->resultTree->SetBranchAddress(paramList[i],&(BB->param[i]));
        BB->resultTree->SetBranchAddress(paramList[i]+"_err",&(BB->paramErr[i]));
        BB->resultTree->GetEntry(0);
      };

      // if MCMC was used, fill param vs sample graphs
      if(minimizer==mkMCMC || minimizer==mkMCMCthenCov) {
        BB->mcmcTree->SetBranchAddress("entry",&entry);
        BB->mcmcTree->SetBranchAddress("nll_MarkovChain_local_",&nll);
        for(int i=0; i<nParams; i++) {
          BB->mcmcTree->SetBranchAddress(paramList[i],&paramval[i]);
        };
        for(int i=0; i<nParams; i++) {
          vTitle = moduList[i] ? moduList[i]->AsymmetryTitle() : "N";
          BB->paramVsSample[i]->SetTitle(
            vTitle+" vs. MCMC sample"/*;sample;"+vTitle*/);
          BB->paramVsSample[i]->GetXaxis()->SetLabelSize(axisTitleSize);
          BB->paramVsSample[i]->GetYaxis()->SetLabelSize(axisTitleSize);
        };
        BB->nllVsSample->SetTitle("-ln(L) vs. MCMC sample");
        BB->nllVsSample->GetXaxis()->SetLabelSize(axisTitleSize);
        BB->nllVsSample->GetYaxis()->SetLabelSize(axisTitleSize);
        for(Long64_t e=0; e<BB->mcmcTree->GetEntries(); e++) {
          BB->mcmcTree->GetEntry(e);
          for(int i=0; i<nParams; i++) {
            BB->paramVsSample[i]->Fill(entry+1,paramval[i]);
          };
          BB->nllVsSample->Fill(entry+1,nll);
        };
      };

    };
  };



  // build graphs and canvases
  TString outfileN;
  TFile * outFile;
  TGraphErrors * paramGr[nParamsMax];
  Int_t cnt;
  TCanvas * paramCanv;
  TCanvas * paramVsSampleCanv;
  TCanvas * cornerCanv;
  TCanvas * autocorrCanv;
  Float_t xMin,xMax,yMin,yMax;
  TLine * zeroLine;
  Int_t nrow,ncol;
  Int_t binListCnt=0;
  TString blStr;
  nextBinList = TObjArrayIter(BruBinSuperList);
  while((BBlist = (TObjArray*) nextBinList())) {

    // setup output file
    blStr = Form("_BL%d",binListCnt); // BL = Bin List
    outfileN = bruDir+Form("/asym_%s%s.root",minimizer_.Data(),blStr.Data());
    outFile = new TFile(outfileN,"RECREATE");


    Tools::PrintSeparator(30,"-");
    printf("%s will have bins:\n",outfileN.Data());
    nextBin = TObjArrayIter(BBlist);
    while((BB = (BruBin*) nextBin())) BB->PrintInfo();
    Tools::PrintSeparator(30,"-");

    // paramter vs. bin mean graphs, for each parameter
    for(int i=0; i<nParams; i++) {

      // define graph
      paramGr[i] = new TGraphErrors();
      paramGr[i]->SetName("gr_"+paramList[i]+blStr);
      vTitle = moduList[i] ? moduList[i]->AsymmetryTitle() : "N";
      paramGr[i]->SetTitle(vTitle+" vs. "+hTitle/*+";"+hTitle+";"+vTitle*/);
      paramGr[i]->GetXaxis()->SetLabelSize(axisTitleSize);
      paramGr[i]->GetYaxis()->SetLabelSize(axisTitleSize);
      paramGr[i]->SetMarkerStyle(kFullCircle);
      paramGr[i]->SetMarkerColor(kAzure);
      paramGr[i]->SetLineColor(kAzure);

      // add points to graph and write
      cnt=0;
      nextBin = TObjArrayIter(BBlist);
      while((BB = (BruBin*) nextBin())) {
        paramGr[i]->SetPoint(cnt,BB->mean,BB->param[i]);
        paramGr[i]->SetPointError(cnt,0,BB->paramErr[i]);
        cnt++;
      };
      paramGr[i]->Write();
    };

    // canvas for paramter vs. bin mean graphs
    ncol=4; nrow=(nParams-1)/ncol+1;
    paramCanv = new TCanvas("canvAsym"+blStr,"canvAsym"+blStr,600*ncol,300*nrow);
    paramCanv->Divide(ncol,nrow);
    for(int i=0; i<nParams; i++) {
      paramCanv->cd(i+1);
      paramCanv->GetPad(i+1)->SetGrid(1,1);
      yMin = asymPlotMin;
      yMax = asymPlotMax;
      if(paramGr[i]->GetYaxis()->GetXmin() < yMin)
        yMin = paramGr[i]->GetYaxis()->GetXmin();
      if(paramGr[i]->GetYaxis()->GetXmax() > yMax)
        yMax = paramGr[i]->GetYaxis()->GetXmax();
      paramGr[i]->GetYaxis()->SetRangeUser(yMin,yMax);
      paramGr[i]->Draw("APE");
      xMin = paramGr[i]->GetXaxis()->GetXmin();
      xMax = paramGr[i]->GetXaxis()->GetXmax();
      if(!paramList[i].Contains("Yld")) {
        zeroLine = new TLine(xMin,0,xMax,0);
        zeroLine->SetLineColor(kBlack);
        zeroLine->SetLineWidth(2);
        zeroLine->SetLineStyle(kDashed);
        zeroLine->Draw();
      };
    };
    paramCanv->Draw();
    paramCanv->Write();

    // parameter vs. sample
    if(minimizer==mkMCMC || minimizer==mkMCMCthenCov) {
      nrow=nParams/ncol+1; // (update for NLL)
      nextBin = TObjArrayIter(BBlist);
      while((BB = (BruBin*) nextBin())) {
        paramVsSampleCanv = new TCanvas(
          Form("paramVsSample_%d",BB->idx)+blStr,
          Form("paramVsSample_%d",BB->idx)+blStr,
          600*ncol,300*nrow);
        paramVsSampleCanv->Divide(ncol,nrow);
        for(int i=0; i<nParams; i++) {
          paramVsSampleCanv->cd(i+1);
          if(logscale) gPad->SetLogx();
          if(!paramList[i].Contains("Yld"))
            BB->paramVsSample[i]->GetYaxis()->SetRangeUser(
              asymPlotMin,asymPlotMax);
          BB->paramVsSample[i]->Draw("HIST");
        };
        paramVsSampleCanv->cd(nParams+1);
        if(logscale) gPad->SetLogx();
        BB->nllVsSample->Draw("HIST");
        paramVsSampleCanv->Write();
        paramVsSampleCanv->Close();

        cornerCanv = (TCanvas*) BB->resultFile->Get("Corner Full Plot")->Clone();
        //cornerCanv = (TCanvas*) BB->resultFile->Get("Corner Plot")->Clone();
        cornerCanv->Write(Form("cornerCanv_%d"+blStr,BB->idx));

        autocorrCanv = (TCanvas*) BB->resultFile->Get("Autocorrelation Plot")->Clone();
        autocorrCanv->Write(Form("autocorrCanv_%d"+blStr,BB->idx));
      };
    };

    binListCnt++;
    outFile->Close();
    printf("produced %s\n",outfileN.Data());
  };

  // cleanup
  nextBinList = TObjArrayIter(BruBinSuperList);
  while((BBlist = (TObjArray*) nextBinList())) {
    nextBin = TObjArrayIter(BBlist);
    while((BB = (BruBin*) nextBin())) {
     BB->resultFile->Close();
    };
  };

};
