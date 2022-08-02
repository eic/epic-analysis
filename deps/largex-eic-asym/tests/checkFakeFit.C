// check fake asymmetry injection
void checkFakeFit(TString infileN="bruspin/asym_minuit_BL0.root") {
  TFile *infile = new TFile(infileN,"READ");
  const int N=5;
  TGraphErrors *gr[N];
  for(int i=0; i<N; i++) gr[i] = (TGraphErrors*)infile->Get(Form("gr_AmpP3I%d_BL0",i));
  Double_t ampVal[N] = { 0.1, 0.1, 0.0, 0.0, 0.0 };
  TF1 *func[N];
  func[0] = new TF1("func0",Form(" x*%f/0.2",ampVal[0]),0,1);
  func[1] = new TF1("func1",Form("-x*%f/0.2",ampVal[1]),0,1);
  func[2] = new TF1("func2","0",0,1);
  func[3] = new TF1("func3","0",0,1);
  func[4] = new TF1("func4","0",0,1);
  TCanvas *canv[N];
  for(int i=0; i<N; i++) {
    canv[i] = new TCanvas(Form("canv%d",i),Form("canv%d",i),600,600);
    canv[i]->SetGrid(1,1);
    canv[i]->SetBottomMargin(0.15);
    canv[i]->SetLeftMargin(0.15);
    func[i]->SetLineColor(kBlack);
    func[i]->SetLineStyle(kDashed);
    gr[i]->Draw("APE");
    gr[i]->GetYaxis()->SetRangeUser(-0.2,0.2);
    gr[i]->Fit("pol1","","",0,1);
    func[i]->Draw("SAME");
  };
};




