// pullWeights: pull the specified weights branch from the tree, and produce
// a `Tweights.root` file, needed to apply weights in brufit.C
//
// - IMPORTANT: execute with `brufit pullWeights.C'(....)'` (you can do `brufit -b -q`, if you want)
//
// - outFileName likely does not have to be "Tweights.root", this is just following
//   the sPlot file name convention in Brufit
//
// 
void pullWeights(
    TString inFileName = "../out/asym.idx.root", // input, indexed tree
    TString weightBranch = "Weight", // name of weight branch
    TString outFileName = "data/Tweights.root", // name of output file
    TString treeName = "tree" // name of TTree in input file
    )
{
  // get input TTree
  TFile *inFile = new TFile(inFileName,"READ");
  TTree *inTree = (TTree*)inFile->Get(treeName);

  // setup Brufit Weights object
  Weights *W = new Weights(weightBranch+"Type");
  W->SetFile(outFileName);
  W->SetSpecies(weightBranch+"Class");
  W->SetIDName("Idx");

  // pull the weights
  W->WeightBySelection(inTree,"",weightBranch);
  W->Save();
  inFile->Close();
};
