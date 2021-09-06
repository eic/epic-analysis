#include "PostProcessor.h"

ClassImp(PostProcessor)

using std::cout;
using std::cerr;
using std::endl;

// constructor
PostProcessor::PostProcessor(
  TString infileN_
)
  : infileN(infileN_)
{

  // settings
  // - summary canvas formatting
  summaryColor[0] = kRed;     summaryStyle[0] = kFullCircle;
  summaryColor[1] = kGreen+1; summaryStyle[1] = kFullTriangleUp;
  summaryColor[2] = kBlue;    summaryStyle[2] = kFullTriangleDown;
  // - general
  gStyle->SetOptStat(0);
  cout << std::fixed; // fixed float precision
  cout << std::setprecision(3);

  // set up input and output files
  cout << "-- postprocess histograms from " << infileN << endl;
  outfileN = infileN;
  outfileN(TRegexp("\\.root$")) = "";
  pngDir = outfileN+".images";
  outfileN += ".canvas.root";
  outfile = new TFile(outfileN,"RECREATE");
  infile = new TFile(infileN,"READ");
  gROOT->ProcessLine(".! mkdir -p "+pngDir);

  // initialize algorithm-specific vars
  this->ResetVars();

};


//=========================================================================
// cleanup and close open files and streams
void PostProcessor::Finish() {
  infile->Close();
  outfile->Close();
  cout << outfileN << " written." << endl;
  cout << pngDir << "/ images written." << endl;
};


//=========================================================================
// reset global variables, such as summary canvases
void PostProcessor::ResetVars() {
  nsum=0;
  ndump=0;
  summaryCanvMap.clear();
  varList.clear();
};


//=========================================================================
/* ALGORITHM: dump histogram counts
 * - output will be appended to `countsFiles`
 * - so far only implemented for 1D
 */
void PostProcessor::DumpHist(TString datFile, TString histSet, TString varName) {
  cout << "dump " << histSet << " : " << varName << " to " << datFile << endl;
  Histos *H = (Histos*) infile->Get(histSet);
  TH1 *hist = H->Hist(varName);
  if(hist->GetDimension()>1) return;
  TString histTformatted = hist->GetTitle();
  histTformatted.ReplaceAll("#in"," in ");
  histTformatted.ReplaceAll("#pm","+-");

  gSystem->RedirectOutput(datFile,"a");
  cout << endl;
  cout << "Histogram: " << histTformatted << endl;

  gSystem->RedirectOutput("tempo","w");
  cout << "bin"
       << " " << varName << "_min"
       << " " << varName << "_max"
       << " " << "counts"
       << " " << "error"
       << " " << endl;
  for(int b=1; b<=hist->GetNbinsX(); b++) {
    cout << b
         << " " << hist->GetBinLowEdge(b)
         << " " << hist->GetBinLowEdge(b+1)
         << " " << hist->GetBinContent(b)
         << " " << hist->GetBinError(b)
         << endl;
  };

  gSystem->RedirectOutput(datFile,"a");
  gROOT->ProcessLine(".! cat tempo | column -t");
  gSystem->RedirectOutput(0);
  gROOT->ProcessLine(".! rm tempo");

};

//=========================================================================

/* ALGORITHM: dump averages of all histograms from `histSet` to the file `datFile`
 * - it is best to call this function in a loop, this function will only dump
 *   one line of information; if it's the first time you called it, it will
 *   dump header information as well
 * - use `cutName` to specify a cut defintion you want to print; e.g., loop
 *   over x bins and set the `cutName` to the x cut, to include columns for x
 *   bin boundaries
 * - when done looping, call `FinishDumpAve(datFile)`
 */
void PostProcessor::DumpAve(TString datFile, TString histSet, TString cutName) {
  cout << "dump averages from " << histSet
       << " to " << datFile << endl;
  Histos *H = (Histos*) infile->Get(histSet);

  // get associated cut, and print header info if it's the first time
  dumpCut = nullptr;
  if(ndump==0) {
    gSystem->RedirectOutput(datFile,"a");
    cout << endl << "Histogram Set:";
  };
  for(CutDef *cut : H->CutDefList) {
    if(cutName.CompareTo(cut->GetVarName(),TString::kIgnoreCase)==0) {
      dumpCut = cut;
    }
    else if(ndump==0) {
      TString cutT = cut->GetCutTitle();
      cutT.ReplaceAll("#in"," in ");
      cutT.ReplaceAll("#pm","+-");
      cout << "   " + cutT;
    };
  };
  if(ndump==0) {
    cout << endl;
    cout << "Averages in bins of " << dumpCut->GetVarTitle() << ":" << endl;
    cout << "------------------------" << endl;
    gSystem->RedirectOutput(0);
  };
  if(dumpCut==nullptr) {
    cerr << "ERROR: in DumpAve, cannot find cut def" << endl;
    return;
  };

  // start temporary file for output; this is so you can use `Columnify` later
  // for pretty print to full `datFile`
  gSystem->RedirectOutput(datFile+".tmp",ndump>0?"a":"w");

  // if this is the first time, populate ordered list of histograms; we keep a
  // local list, in case you call this function on a `Histos` with different
  // ordering; we want to maintain order for printout
  if(ndump==0) {
    // header for bin columns
    cout << "bin"
         << " " << dumpCut->GetVarTitle() << "_min"
         << " " << dumpCut->GetVarTitle() << "_max"
         << " counts";
    // loop over 1D histograms, build `varList`,
    // and print header for ave columns
    for(TString varName : H->VarNameList) {
      if(H->Hist(varName)->GetDimension()==1) {
        if(varName.Contains("xsec")) continue;
        varList.push_back(varName);
        cout << " <" + varName + ">";
      };
    };
    // header newline
    cout << endl;
  };

  // loop over local list of variables, and print their averages
  TH1 *hist;
  cout << ndump
       << " " << dumpCut->GetMin()
       << " " << dumpCut->GetMax();
  Bool_t first=true;
  Double_t counts;
  for(TString varName : varList) {
    hist = H->Hist(varName);
    if(first) {
      counts = hist->GetEntries();
      cout << " " << counts;
      first=false;
    }
    else if(hist->GetEntries()!=counts) {
      cerr << "WARNING: mismatch counts" << endl;
    };
    cout << " " << hist->GetMean();
  };
  cout << endl;
  gSystem->RedirectOutput(0);


  // iterate row counter
  ndump++;

};
void PostProcessor::FinishDumpAve(TString datFile) {
  this->Columnify(datFile+".tmp",datFile);
  gROOT->ProcessLine(".! rm "+datFile+".tmp");
  this->ResetVars();
};

//=========================================================================

/* ALGORITHM: draw a single histogram to a canvas, and write it
 * - since `histSet` names can be hard to read, you can use `outName` to give a
 *   "nickname" to `histSet`, which the canvas name will include
 */
void PostProcessor::DrawSingle(TString histSet, TString histName) {

  Histos *H = (Histos*) infile->Get(histSet);
  TH1 *hist = H->Hist(histName);

  TString canvN = "canv_"+histName+"_"+H->GetSetName();
  TCanvas *canv = new TCanvas(canvN,canvN, dimx, dimy);

  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(kFullCircle);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  
  // determine draw type (TODO: probably could be generalized somehow)
  TString drawStr = "";
  switch(hist->GetDimension()) {
    case 1:
      drawStr = "EX0 P";
      if(histName=="Q_xsec") {
        hist->GetXaxis()->SetRangeUser(1,10);
        hist->GetYaxis()->SetRangeUser(1e-6,1);
        drawStr = "E P";
      };
      break;
    case 2:
      drawStr = "COLZ";
      break;
    case 3:
      drawStr = "BOX";
      break;
  };

  hist->Draw(drawStr);

  canv->SetGrid(1,1);
  canv->SetLogx(H->GetHistConfig(histName)->logx);
  canv->SetLogy(H->GetHistConfig(histName)->logy);
  canv->SetLogz(H->GetHistConfig(histName)->logz);
  canv->SetBottomMargin(0.15);
  canv->SetLeftMargin(0.15);
  canv->Print(pngDir+"/"+canvN+".png");
  outfile->cd("/");
  canv->Write();
  outfile->cd("/");

  /* // deprecated, for combining single plots; TODO if needed, make a separate method
  TH1 *histClone = (TH1*) hist->Clone();
  histClone->SetLineColor  (nsum<nsumMax?summaryColor[nsum]:kBlack);
  histClone->SetMarkerColor(nsum<nsumMax?summaryColor[nsum]:kBlack);
  histClone->SetMarkerStyle(nsum<nsumMax?summaryStyle[nsum]:kFullCircle);

  if(nsum==0) {
    summaryCanv = new TCanvas(
        "summaryCanv_"+histName,
        "summaryCanv_"+histName,
        dimx,dimy
        );
    summaryCanv->SetGrid(1,1);
    summaryCanv->SetLogx(H->GetHistConfig(histName)->logx);
    summaryCanv->SetLogy(H->GetHistConfig(histName)->logy);
    summaryCanv->SetLogz(H->GetHistConfig(histName)->logz);
    summaryCanv->SetBottomMargin(0.15);
    summaryCanv->SetLeftMargin(0.15);
  };
  summaryCanv->cd();
  histClone->Draw(drawStr+(nsum>0?" SAME":""));
  nsum++;
  */
};


//=========================================================================

/* ALGORITHM: draw a ratio of all 1D histograms in the specified histogram set
* - the ratio will be of `numerSet` over `denomSet`
* - canvases will be created, and depending on the setting of `plotRatioOnly`,
*   either just the ratio will be plotted, or the two histograms will also be
*   plotted
* - Use `outName` to specify a names for output canvases
* - summary canvases are also created, which can combine multiple ratio plots;
*   they are accumulated in `summaryCanvMap` and colors are set by `nsum`; call
*   `ResetVars` to reset the accumulation
 * - when done looping, call `FinishDrawRatios`
*/
void PostProcessor::DrawRatios(
    TString outName, TString numerSet, TString denomSet, Bool_t plotRatioOnly
) {

  cout << "draw ratios " << outName << "..." << endl;
  enum HHenum {num,den};
  Histos *HH[2];
  HH[num] = (Histos*) infile->Get(numerSet);
  HH[den] = (Histos*) infile->Get(denomSet);
  TH1 *hist[2];
  TH1D *ratio;
  Double_t ny,ry,err;
  TCanvas *canv;
  TString histT[2], tok[2], uniqT[2], sameT[2];
  Ssiz_t tf;

  outfile->cd("/");
  outfile->mkdir(outName);
  outfile->cd(outName);

  // loop over 1D histograms
  for(TString varName : HH[num]->VarNameList) {
    hist[num] = HH[num]->Hist(varName);
    if(hist[num]->GetDimension()==1) {
      hist[den] = HH[den]->Hist(varName);

      // filter title
      for(int f=0; f<2; f++)  histT[f] = hist[f]->GetTitle();
      for(int f=0; f<2; f++) {
        tf = 0;
        sameT[f] = "";
        uniqT[f] = "";
        while(histT[f].Tokenize(tok[f],tf,", ")) {
          if(histT[(f+1)%2].Contains(tok[f])) sameT[f] += ", " + tok[f];
          else                                uniqT[f] += ", " + tok[f];
        };
        sameT[f](TRegexp("^, "))="";
        uniqT[f](TRegexp("^, "))="";
        /*cout <<  "histT = " << histT[f] << endl;
          cout <<  "uniqT = " << uniqT[f] << endl;
          cout <<  "sameT = " << sameT[f] << endl;*/
        hist[f]->SetTitle(sameT[f]);
      };
      TString ratioStr = "ratio("+uniqT[num]+"/"+uniqT[den]+")";

      // calculate ratio
      ratio = (TH1D*) hist[num]->Clone();
      ratio->Divide(hist[den]);
      TString ratioT = ratio->GetTitle();
      TString ratioSummaryT = ratioT;
      ratioT(TRegexp("distribution")) = ratioStr;
      ratioSummaryT(TRegexp("distribution")) = "ratios";
      ratio->SetTitle(ratioT);

      // calculate ratio errors, assuming ratio represents an efficiency
      for(int b=0; b<=ratio->GetNbinsX(); b++) {
        ny = hist[den]->GetBinContent(b);
        ry = ratio->GetBinContent(b);
        err = TMath::Abs(ny)>0 ?
          TMath::Sqrt( ry*(1-ry) / ny )
          : 0;
        ratio->SetBinError(b,err);
      };

      // canvas
      TString canvN = "canv_"+outName+"_"+varName;
      Int_t dimx=800;
      Int_t dimy=700;
      canv = new TCanvas(canvN,canvN, dimx*(plotRatioOnly?1:2), dimy);

      // draw plots
      int numPads;
      if(!plotRatioOnly) { canv->Divide(2,1); numPads=2; }
      else { canv->Divide(1,1); numPads=1; };
      if(!plotRatioOnly) {
        canv->cd(2);
        hist[num]->SetLineColor(kYellow-8);
        hist[den]->SetLineColor(kAzure-7);
        hist[num]->SetMarkerColor(kYellow-8);
        hist[den]->SetMarkerColor(kAzure-7);
        hist[num]->SetMarkerStyle(kFullCircle);
        hist[den]->SetMarkerStyle(kFullCircle);
        for(int f=0; f<2; f++) {
          hist[f]->SetMarkerSize(1.0);
          hist[f]->SetLineWidth(2);
          hist[f]->GetXaxis()->SetLabelSize(0.06);
          hist[f]->GetYaxis()->SetLabelSize(0.06);
          hist[f]->GetXaxis()->SetTitleSize(0.06);
          hist[f]->GetYaxis()->SetTitleSize(0.06);
          hist[f]->GetXaxis()->SetTitleOffset(1.2);
          //if(plotName=="PperpDistLin") hist[f]->GetXaxis()->SetRangeUser(0,1.7);
        };
        hist[den]->Draw("EX0 P");
        hist[num]->Draw("EX0 P SAME");
      };
      canv->cd(1);
      ratio->GetYaxis()->SetTitle(ratioStr);
      ratio->SetLineColor(kBlack);
      ratio->SetMarkerColor(kBlack);
      ratio->SetMarkerStyle(kFullCircle);
      ratio->SetMarkerSize(1.0);
      ratio->SetLineWidth(2);
      ratio->GetXaxis()->SetLabelSize(0.06);
      ratio->GetYaxis()->SetLabelSize(0.06);
      ratio->GetXaxis()->SetTitleSize(0.06);
      ratio->GetYaxis()->SetTitleSize(0.06);
      ratio->GetXaxis()->SetTitleOffset(1.2);
      ratio->GetYaxis()->SetRangeUser(0,1.3);
      ratio->Draw("EX0 P");
      for(int pad=1; pad<=numPads; pad++) {
        canv->GetPad(pad)->SetGrid(1,1);
        canv->GetPad(pad)->SetLogx(HH[num]->GetHistConfig(varName)->logx);
        if(pad>1) {
          canv->GetPad(pad)->SetLogy(HH[num]->GetHistConfig(varName)->logy);
          canv->GetPad(pad)->SetLogz(HH[num]->GetHistConfig(varName)->logz);
        };
        canv->GetPad(pad)->SetBottomMargin(0.15);
        canv->GetPad(pad)->SetLeftMargin(0.15);
      };
      canv->Print(pngDir+"/"+canvN+".png");
      canv->Write();

      // add to summary canvas
      TH1D *ratioSummary = (TH1D*) ratio->Clone();
      ratioSummary->SetTitle(ratioSummaryT);
      ratioSummary->SetYTitle("ratio");
      ratioSummary->SetLineColor  (nsum<nsumMax?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerColor(nsum<nsumMax?summaryColor[nsum]:kBlack);
      ratioSummary->SetMarkerStyle(nsum<nsumMax?summaryStyle[nsum]:kFullCircle);
      auto kv = summaryCanvMap.find(varName);
      if(kv==summaryCanvMap.end()) {
        summaryCanv = new TCanvas(
            "summaryCanv_ratio"+varName,
            "summaryCanv_ratio"+varName,
            dimx,dimy
            );
        summaryCanv->SetGrid(1,1);
        summaryCanv->SetLogx(HH[num]->GetHistConfig(varName)->logx);
        summaryCanv->SetBottomMargin(0.15);
        summaryCanv->SetLeftMargin(0.15);
        summaryCanvMap.insert(std::pair<TString,TCanvas*>(varName,summaryCanv));
        ratioSummary->Draw();
      } else {
        kv->second->cd();
        ratioSummary->Draw("SAME");
      };

    };
  };
  outfile->cd("/");
  nsum++;
};
void PostProcessor::FinishDrawRatios(TString summaryDir) {
  // write summary canvas
  outfile->cd("/");
  outfile->mkdir(summaryDir);
  outfile->cd(summaryDir);
  for(auto const& kv : summaryCanvMap) {
    kv.second->Write();
    kv.second->Print(pngDir+"/"+summaryDir+"_"+kv.first+".png");
  };
  outfile->cd("/");
  this->ResetVars();
};



//=========================================================================
// TEXT FILE: start new text file, append to text file, print text file
void PostProcessor::StartTextFile(TString datFile, TString firstLine) {
  gSystem->RedirectOutput(datFile,"w");
  if(firstLine!="") cout << firstLine << endl;
  gSystem->RedirectOutput(0);
};
void PostProcessor::AppendToTextFile(TString datFile, TString appendText="") {
  gSystem->RedirectOutput(datFile,"a");
  cout << appendText << endl;
  gSystem->RedirectOutput(0);
};
void PostProcessor::PrintTextFile(TString datFile) {
  gROOT->ProcessLine(".! cat "+datFile);
};


//=========================================================================
// TEXT FILE: pipe `inputFile` through `column -t` and append output to `outputFile`
void PostProcessor::Columnify(TString inputFile, TString outputFile) {
  gSystem->RedirectOutput(outputFile,"a");
  gROOT->ProcessLine(".! cat "+inputFile+" | column -t");
  gSystem->RedirectOutput(0);
};


//=========================================================================
// get pointer to BinSet associated with variable
BinSet *PostProcessor::GetBinSet(TString varName) {
  return (BinSet*) infile->Get(varName+"_bins");
};
// get pointer to CutDef associated with a particular bin
CutDef *PostProcessor::GetBinCut(TString varName, Int_t binNum) {
  return this->GetBinSet(varName)->Cut(binNum);
};
// get list of bin numbers
std::vector<int> PostProcessor::GetBinNums(TString varName) {
  BinSet *B = this->GetBinSet(varName);
  std::vector<int> retVec;
  if(B==nullptr) {
    cerr << "ERROR: unknown variable " << varName << " in PostProcessor::GetBinNums" << endl;
  } else {
    for(int n=0; n<B->GetBinList()->GetEntries(); n++) retVec.push_back(n);
  };
  return retVec;
};

// return true if the bin is "full" range, and it's not the only bin
Bool_t PostProcessor::SkipFull(TString varName, Int_t binNum) {
  return GetBinSet(varName)->GetNumBins() > 1 &&
         GetBinCut(varName,binNum)->GetCutType()=="Full";
};


//=========================================================================
// destructor
PostProcessor::~PostProcessor() {
};

