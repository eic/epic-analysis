#include "Analysis.h"

ClassImp(Analysis)

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
Analysis::Analysis(
  TString infile_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_
)
  : infile(infile_)
  , eleBeamEn(eleBeamEn_)
  , ionBeamEn(ionBeamEn_)
  , crossingAngle(crossingAngle_)
{
  // set bin schemes
  AddBinScheme("pt","p_{T}");
  AddBinScheme("z","z");
  AddBinScheme("x","x");
  AddBinScheme("q","Q");
  AddBinScheme("y","y");
  // set final states
  // TODO: consider making this another bin scheme
  PIDtoEnum.insert(std::pair<int,int>(211,pPip));
  //PIDtoEnum.insert(std::pair<int,int>(-211,pPim));
  // TODO: generalized diagonalizer
  diagonalPtXZ = false;
  diagonalXZQ = false;
};


//=============================================
// perform the analysis
//=============================================
void Analysis::Execute() {

  // read delphes tree
  cout << "-- running analysis of " << infile << endl;
  TChain *chain = new TChain("Delphes");
  chain->Add(infile);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  Long64_t ENT = tr->GetEntries();

  // define output file
  TString outfileN = infile;
  outfileN(TRegexp("^.*/")) = "";
  outfileN = "out/histos."+outfileN;
  TFile *outfile = new TFile(outfileN,"RECREATE");

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itPIDSystemsTrack(tr->UseBranch("PIDSystemsTrack"));

  // number of bins
  const Int_t NptBins = BinScheme("pt")->GetBinList()->GetEntries();
  const Int_t NxBins = BinScheme("x")->GetBinList()->GetEntries();
  const Int_t NzBins = BinScheme("z")->GetBinList()->GetEntries();
  const Int_t NqBins = BinScheme("q")->GetBinList()->GetEntries();
  const Int_t NyBins = BinScheme("y")->GetBinList()->GetEntries();

  // sets of histogram sets
  // - `histSet` is a data structure for storing and organizing pointers to
  //   sets of histograms (`Histos` objects)
  // - `histSet*List` are used as temporary lists of relevant `Histos` pointers
  // - TODO: if we add one more dimension, 7D array will probably break; need
  //         better data structure
  Histos *histSet[NptBins][NxBins][NzBins][NqBins][NyBins][NPart];
  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetFillList;
  std::vector<int> v_pt, v_x, v_z, v_q, v_y;
  // instantiate Histos sets, and populate 
  TString plotN,plotT;
  cout << "Define histograms..." << endl;
  for(int bpt=0; bpt<NptBins; bpt++) { // - loop over pT bins
    for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
      for(int bz=0; bz<NzBins; bz++) { // - loop over z bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q bins
          if(CheckDiagonal(bpt,bx,bz,bq)) continue;
          for(int by=0; by<NyBins; by++) { // - loop over y bins

            // set plot name
            plotN  = "_" + BinScheme("pt")->Cut(bpt)->GetVarName() + Form("%d",bpt);
            plotN += "_" + BinScheme("x")->Cut(bx)->GetVarName() + Form("%d",bx);
            plotN += "_" + BinScheme("z")->Cut(bz)->GetVarName() + Form("%d",bz);
            plotN += "_" + BinScheme("q")->Cut(bq)->GetVarName() + Form("%d",bq);
            plotN += "_" + BinScheme("y")->Cut(by)->GetVarName() + Form("%d",by);

            // set plot title
            plotT  = ", " + BinScheme("pt")->Cut(bpt)->GetCutTitle();
            plotT += ", " + BinScheme("x")->Cut(bx)->GetCutTitle();
            plotT += ", " + BinScheme("z")->Cut(bz)->GetCutTitle();
            plotT += ", " + BinScheme("q")->Cut(bq)->GetCutTitle();
            plotT += ", " + BinScheme("y")->Cut(by)->GetCutTitle();

            // loop over particles
            histSet[bpt][bx][bz][bq][by][pPip] = new Histos("pipTrack"+plotN,"#pi^{+} tracks"+plotT);
            //histSet[bpt][bx][bz][bq][by][pPim] = new Histos("pimTrack"+plotN,"#pi^{-} tracks"+plotT);

            // define set of histograms for this bin
            for(int bp=0; bp<NPart; bp++) {
              HS = histSet[bpt][bx][bz][bq][by][bp];

              // HISTOGRAMS ================================================
              // -- DIS kinematics
              HS->DefineHist2D("Q2vsX","x","Q^{2}","","GeV^{2}",
                  NBINS,1e-3,1,
                  NBINS,1,100,
                  true,true
                  );
              HS->DefineHist1D("Q","Q","GeV",NBINS,1.0,11.0,true,true);
              HS->DefineHist1D("x","x","",NBINS,1e-3,1.0,true,true);
              HS->DefineHist1D("y","y","",NBINS,1e-5,1,true);
              HS->DefineHist1D("W","W","GeV",NBINS,0,15);
              // -- hadron 4-momentum
              HS->DefineHist1D("pLab","p_{lab}","GeV",NBINS,0,10);
              HS->DefineHist1D("pTlab","p_{T}^{lab}","GeV",NBINS,1e-2,3,true);
              HS->DefineHist1D("etaLab","#eta_{lab}","",NBINS,-5,5);
              HS->DefineHist1D("phiLab","#phi_{lab}","",NBINS,-TMath::Pi(),TMath::Pi());
              // -- hadron kinematics
              HS->DefineHist1D("z","z","",NBINS,0,1);
              HS->DefineHist1D("pT","p_{T}","GeV",NBINS,1e-2,3,true);
              HS->DefineHist1D("qT","q_{T}","GeV",NBINS,1e-2,5,true);
              HS->DefineHist1D("qTq","q_{T}/Q","",NBINS,1e-2,3,true);
              HS->DefineHist1D("mX","m_{X}","GeV",NBINS,0,20);
              HS->DefineHist1D("phiH","#phi_{h}","",NBINS,-TMath::Pi(),TMath::Pi());
              HS->DefineHist1D("phiS","#phi_{S}","",NBINS,-TMath::Pi(),TMath::Pi());
              // -- cross sections
              //HS->DefineHist1D("Q_xsec","Q","GeV",10,0.5,10.5,false,true); // linear
              HS->DefineHist1D("Q_xsec","Q","GeV",10,1.0,10.0,true,true); // log
              HS->Hist("Q_xsec")->SetMinimum(1e-10);
              // ===========================================================

              // store cut definitions with histogram sets, then add histogram sets full list
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef(BinScheme("pt")->Cut(bpt));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef(BinScheme("x")->Cut(bx));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef(BinScheme("z")->Cut(bz));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef(BinScheme("q")->Cut(bq));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef(BinScheme("y")->Cut(by));
              histSetList.push_back(histSet[bpt][bx][bz][bq][by][bp]);
            };

          };
        };
      };
    };
  };


  // calculate integrated luminosity
  // - cross sections are hard-coded, coped from pythia output
  Int_t eleBeamEnInt = (Int_t) eleBeamEn;
  Int_t ionBeamEnInt = (Int_t) ionBeamEn;
  Double_t xsecTot; // [nb]
  if     (eleBeamEnInt==5  && ionBeamEnInt==41 ) xsecTot=297.9259;
  else if(eleBeamEnInt==18 && ionBeamEnInt==275) xsecTot=700.0; // TODO; this is approximate
  else {
    cerr << "WARNING: unknown cross section; integrated lumi will be wrong" << endl;
    xsecTot=1;
  };
  Long64_t numGen = tr->GetEntries();
  Double_t lumi = numGen/xsecTot; // [nb^-1]
  TString sep = "--------------------------------------------";
  cout << sep << endl;
  cout << "assumed total cross section: " << xsecTot << " nb" << endl;
  cout << "number of generated events:  " << numGen << endl;
  cout << "Integrated Luminosity:       " << lumi << "/nb" << endl;
  cout << sep << endl;


  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);


  // vars
  Double_t eleP,maxEleP;
  int pid,bpart;


  // event loop =========================================================
  //ENT = 1000; // limiter
  cout << "begin event loop..." << endl;
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%100000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
    tr->ReadEntry(e);

    // electron loop
    // - finds max-momentum electron
    itElectron.Reset();
    maxEleP = 0;
    while(Electron *ele = (Electron*) itElectron()) {
      eleP = ele->PT * TMath::CosH(ele->Eta);
      if(eleP>maxEleP) {
        maxEleP = eleP;
        kin->vecElectron.SetPtEtaPhiM(
            ele->PT,
            ele->Eta,
            ele->Phi,
            Kinematics::ElectronMass()
            );
      };
    };
    if(maxEleP<0.001) continue; // no scattered electron found

    // get hadronic final state variables
    kin->GetHadronicFinalState(itTrack, itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itPIDSystemsTrack, itParticle);

    // calculate DIS kinematics
    kin->CalculateDISbyElectron();



    // track loop
    itTrack.Reset();
    while(Track *trk = (Track*) itTrack()) {
      //cout << e << " " << trk->PID << endl;

      // - check PID, to see if it's a particle we're interested in for
      //   histograms; if not, proceed to next
      pid = trk->PID;
      auto kv = PIDtoEnum.find(pid);
      if(kv!=PIDtoEnum.end()) bpart = kv->second;
      else continue;


      // get parent particle, to check if pion is from vector meson
      GenParticle *trkParticle = (GenParticle*)trk->Particle.GetObject();
      TObjArray *brParticle = (TObjArray*)itParticle.GetCollection();
      GenParticle *parentParticle = (GenParticle*)brParticle->At(trkParticle->M1);
      int parentPID = (parentParticle->PID);


      // calculate hadron kinematics
      kin->vecHadron.SetPtEtaPhiM(
          trk->PT,
          trk->Eta,
          trk->Phi,
          trk->Mass /* TODO: do we use track mass here ?? */
          );
      kin->CalculateHadronKinematics();

      // apply cuts
      if(kin->CutFull()) {

        // decide which histogram sets to fill
        // -- `v_A` will be the list of bins that this event's variable `A`
        //    will be a part of; we use these lists to determine the list
        //    of histogram sets to fill
        // - check pT bin
        CheckBins( BinScheme("pt"), v_pt, kin->pT );
        CheckBins( BinScheme("x"),  v_x,  kin->x );
        CheckBins( BinScheme("z"),  v_z,  kin->z );
        CheckBins( BinScheme("q"),  v_q,  TMath::Sqrt(kin->Q2) );
        CheckBins( BinScheme("y"),  v_y,  kin->y );


        // build list of histogram sets to fill
        histSetFillList.clear();
        for(int bpt : v_pt) {
          for(int bx : v_x) {
            for(int bz : v_z) {
              for(int bq : v_q) {
                for(int by : v_y) {
                  if(!CheckDiagonal(bpt,bx,bz,bq)) {
                    histSetFillList.push_back(histSet[bpt][bx][bz][bq][by][bpart]);
                  };
                };
              };
            };
          };
        };

        // loop through list of histogram sets, and fill them
        for(Histos *H : histSetFillList) {
          // DIS kinematics
          H->Hist("Q2vsX")->Fill(kin->x,kin->Q2);
          H->Hist("Q")->Fill(TMath::Sqrt(kin->Q2));
          H->Hist("x")->Fill(kin->x);
          H->Hist("W")->Fill(kin->W);
          H->Hist("y")->Fill(kin->y);
          // hadron 4-momentum
          H->Hist("pLab")->Fill(kin->pLab);
          H->Hist("pTlab")->Fill(kin->pTlab);
          H->Hist("etaLab")->Fill(kin->etaLab);
          H->Hist("phiLab")->Fill(kin->phiLab);
          // hadron kinematics
          H->Hist("z")->Fill(kin->z);
          H->Hist("pT")->Fill(kin->pT);
          H->Hist("qT")->Fill(kin->qT);
          H->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2));
          H->Hist("mX")->Fill(kin->mX);
          H->Hist("phiH")->Fill(kin->phiH);
          H->Hist("phiS")->Fill(kin->phiS);
          // cross sections
          H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),1.0/lumi);
        };

      };
    };
  };
  cout << "end event loop" << endl;
  // event loop end =========================================================



  // print yields in each bin
  cout << sep << endl << "Histogram Entries:" << endl;
  for(Histos *H : histSetList) {
    cout << H->GetSetTitle() << " ::: "
         << H->Hist("Q2vsX")->GetEntries()
         << endl;
  };



  // write histograms
  cout << sep << endl;
  outfile->cd();
  for(Histos *H : histSetList) H->WriteHists(outfile);
  for(Histos *H : histSetList) H->Write();

  // write binning schemes
  for(auto const &kv : binSchemes) kv.second->Write(kv.first+"_bins");

  // close output
  outfile->Close();
  cout << outfileN << " written." << endl;

  // call draw program
  /*
  TString cmd = "./draw.exe "+outfileN;
  system(cmd.Data());
  */
};



//=============================================



// access bin scheme by name
BinSet *Analysis::BinScheme(TString varname) {
  BinSet *ret;
  try { ret = binSchemes.at(varname); }
  catch(const std::out_of_range &ex) {
    cerr << "ERROR: bin scheme "
         << varname << " not found" << endl;
    return nullptr;
  };
  return ret;
};

// add a new bin scheme
void Analysis::AddBinScheme(TString varname, TString vartitle) {
  binSchemes.insert(
    std::pair<TString,BinSet*>(varname,new BinSet(varname,vartitle))
    );
  // TODO: for now, we need to have at least one bin in each dimension,
  // otherwise for loops won't run; when we generalize the `histSet` data
  // structure, hopefully we can also drop this requirement; the current
  // workaround is to add a `full` bin to each dimension
  BinScheme(varname)->BuildBin("Full");
};
 
// return true, if a diagonal mode is on and this is an 
// off-diagonal bin; if a diagonal mode is not on, always 
// return true
Bool_t Analysis::CheckDiagonal(int cpt, int cx, int cz, int cq) {
  if(diagonalPtXZ) return ( cpt!=cx || cx!=cz );
  else if(diagonalXZQ) return ( cx!=cz || cz!=cq );
  else return true;
};

// scan through bin set `bs`, checking each one; the vector `v` will
// contain the list of array indices for which the cut on `var` is satisfied
void Analysis::CheckBins(BinSet *bs, std::vector<int> &v, Double_t var) {
  v.clear();
  for(int b=0; b<bs->GetBinList()->GetEntries(); b++) {
    if(bs->Cut(b)->CheckCut(var)) v.push_back(b);
  };
};

// destructor
Analysis::~Analysis() {
};

