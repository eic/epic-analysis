// computes quantiles for a variety of distributions in plots.root
// -- use `constraint` string to specify a constraint, e.g., `abs(Mh-0.77)<0.1`;
//    this is useful for constructing multidimensional bins

const Int_t N = 100;
Double_t q[N];
Double_t p[N];
int i;
Int_t Nquant;
TString constraint;

void PrintNums(TString name);
TFile * f;
TTree * t;
TH1D * d;
TLine * l;

void quantiles(
  Int_t Nquant_=6,
  TString constraint_="",
  TString treeFile="../out/test.simple.tree.example_5x41.root"
) {
  Nquant = Nquant_;
  constraint = constraint_;
  if(Nquant>N) { fprintf(stderr,"ERROR: Nquant too big\n"); return; }
  f = new TFile(treeFile,"READ");
  t = (TTree*) f->Get("tree");

  PrintNums("X");
  PrintNums("Q2");
  PrintNums("Z");
  PrintNums("PhPerp");
};


void PrintNums(TString name) {

  Double_t min = t->GetMinimum(name);
  Double_t max = t->GetMaximum(name);
  Double_t minTmp = min;
  min -= abs(max-min)*0.05;
  max += abs(max-min)*0.05;

  d = new TH1D(name+"_dist",name+" quantiles",1000,min,max);
  t->Project(name+"_dist",name,constraint);
  new TCanvas();
  d->Draw();

  for(i=0; i<Nquant; i++) p[i] = Double_t(i+1)/Nquant;
  d->GetQuantiles(Nquant,q,p);

  for(i=0; i<Nquant-1; i++) {
    l = new TLine(q[i],0,q[i],d->GetMaximum());
    l->SetLineWidth(3);
    l->Draw();
  };

  // print code for brufit.C
  printf("\n");
  printf("Double_t %sbins[%d] = { ",name.Data(),Nquant+1);
  printf("%.3f",minTmp);
  for(i=0; i<Nquant; i++) printf(", %.3f",q[i]);
  printf(" };\n");
  printf("B->Bin(\"%s\",%d,%sbins);\n",name.Data(),Nquant,name.Data());
  printf("\n");
};
