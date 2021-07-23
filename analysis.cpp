#include <stdlib.h>

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

// subroutines
Bool_t CheckDiagonal(int cpt, int cx, int cz);
void CheckBins(TObjArray *binArr, std::vector<int> &v, Double_t var);

// globals
Bool_t diagonalBinsOnly;

//=========================================================================

int main(int argc, char **argv) {

  // ARGUMENTS ////////////////////////////////////////////////
  TString infile;
  Double_t eleBeamEn = 5; // GeV
  Double_t ionBeamEn = 41; // GeV
  Double_t crossingAngle = 0; // mrad
  if(argc<=1) {
    cout << "USAGE: " << argv[0]
         << " [rootfile with Delphes tree]"
         << " [eleBeamEn(def=" << eleBeamEn << " GeV)]"
         << " [ionBeamEn(def=" << ionBeamEn << " GeV)]"
         << " [crossingAngle(def=" << crossingAngle << " mrad)]"
         << endl;
         return 1;
  };
  if(argc>1) infile = TString(argv[1]);
  if(argc>2) eleBeamEn = (Double_t)strtof(argv[2],NULL);
  if(argc>3) ionBeamEn = (Double_t)strtof(argv[3],NULL);
  if(argc>4) crossingAngle = (Double_t)strtof(argv[4],NULL);
  

  // =========================================
  // BINNING
  // =========================================

  // arrays of bins
  TObjArray *ptBins = new TObjArray();
  TObjArray *xBins = new TObjArray();
  TObjArray *zBins = new TObjArray();
  TObjArray *qBins = new TObjArray();
  TObjArray *yBins = new TObjArray();

  // if this is true, only take 'diagonal' elements of the multi
  // dimensional array of possible pT,x,z bins; this is useful
  // if you want to check specific bins
  diagonalBinsOnly = false; // default value, may be overriden below

  // full bins; useful to have these first -----------------
  ptBins->AddLast(new CutDef("pt","p_{T}","Full"));
  xBins->AddLast(new CutDef("x","x","Full"));
  zBins->AddLast(new CutDef("z","z","Full"));
  qBins->AddLast(new CutDef("q","Q","Full"));
  yBins->AddLast(new CutDef("y","y","Full"));


  // cross check cross section -------------------
  ///*
  diagonalBinsOnly = true;
  // slide 11
  xBins->AddLast(new CutDef("x","x","CenterDelta", 0.3, 0.05 ));
  zBins->AddLast(new CutDef("z","z","CenterDelta", 0.7, 0.05 ));
  ptBins->AddLast(new CutDef("pt","p_{T}","CenterDelta", 0.5, 0.05 ));
  // slide 13
  xBins->AddLast(new CutDef("x","x","CenterDelta", 0.6, 0.05 ));
  zBins->AddLast(new CutDef("z","z","CenterDelta", 0.5, 0.05 ));
  ptBins->AddLast(new CutDef("pt","p_{T}","CenterDelta", 0.55, 0.05 ));
  // slide 14
  xBins->AddLast(new CutDef("x","x","CenterDelta", 0.1, 0.05 ));
  zBins->AddLast(new CutDef("z","z","CenterDelta", 0.7, 0.05 ));
  ptBins->AddLast(new CutDef("pt","p_{T}","CenterDelta", 0.15, 0.05 ));
  //*/

  // cross check dihadrons from EIC smear ------------
  /*
  zBins->AddLast(new CutDef("z","z","Range", 0.2, 0.3 ));
  zBins->AddLast(new CutDef("z","z","Range", 0.3, 0.9 ));
  */

  // - Q bins ------------------------------------
  ///*
  BinSet *qBinScheme = new BinSet("q","Q",10,1,11,false);
  qBins = qBinScheme->bins; // (overwrites any previous qBins)
  //*/

  // - y-minimum cuts ----------------------------
  ///*
  //yBins->AddLast(new CutDef("y","y","Min",0.03));
  yBins->AddLast(new CutDef("y","y","Min",0.05));
  //yBins->AddLast(new CutDef("y","y","Min",0.10));
  //*/

  // - particle species --------------------------
  std::map<int,int> PIDtoEnum;
  enum partEnum{
    pPip,
    //pPim,
    NPart
  };
  PIDtoEnum.insert(std::pair<int,int>(211,pPip));
  //PIDtoEnum.insert(std::pair<int,int>(-211,pPim));


  // ============================================
  const Int_t NptBins = ptBins->GetEntries();
  const Int_t NxBins = xBins->GetEntries();
  const Int_t NzBins = zBins->GetEntries();
  const Int_t NqBins = qBins->GetEntries();
  const Int_t NyBins = yBins->GetEntries();

  /////////////////////////////////////////////////////////////



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
          if(CheckDiagonal(bpt,bx,bz)) continue;
          for(int by=0; by<NyBins; by++) { // - loop over y bins

            // set plot name
            plotN  = "_" + ((CutDef*)ptBins->At(bpt))->GetVarName() + Form("%d",bpt);
            plotN += "_" + ((CutDef*)xBins->At(bx))->GetVarName() + Form("%d",bx);
            plotN += "_" + ((CutDef*)zBins->At(bz))->GetVarName() + Form("%d",bz);
            plotN += "_" + ((CutDef*)qBins->At(bq))->GetVarName() + Form("%d",bq);
            plotN += "_" + ((CutDef*)yBins->At(by))->GetVarName() + Form("%d",by);

            // set plot title
            plotT  = ", " + ((CutDef*)ptBins->At(bpt))->GetCutTitle();
            plotT += ", " + ((CutDef*)xBins->At(bx))->GetCutTitle();
            plotT += ", " + ((CutDef*)zBins->At(bz))->GetCutTitle();
            plotT += ", " + ((CutDef*)qBins->At(bq))->GetCutTitle();
            plotT += ", " + ((CutDef*)yBins->At(by))->GetCutTitle();

            // loop over particles
            histSet[bpt][bx][bz][bq][by][pPip] = new Histos("pipTrack"+plotN,"#pi^{+} tracks"+plotT);
            //histSet[bpt][bx][bz][bq][by][pPim] = new Histos("pimTrack"+plotN,"#pi^{-} tracks"+plotT);

            // store cut definitions with histogram sets, then add histogram sets full list
            for(int bp=0; bp<NPart; bp++) {
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef((CutDef*)ptBins->At(bpt));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef((CutDef*)xBins->At(bx));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef((CutDef*)zBins->At(bz));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef((CutDef*)qBins->At(bq));
              histSet[bpt][bx][bz][bq][by][bp]->AddCutDef((CutDef*)yBins->At(by));
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
        CheckBins( ptBins, v_pt, kin->pT );
        CheckBins( xBins,  v_x,  kin->x );
        CheckBins( zBins,  v_z,  kin->z );
        CheckBins( qBins,  v_q,  TMath::Sqrt(kin->Q2) );
        CheckBins( yBins,  v_y,  kin->y );


        // build list of histogram sets to fill
        histSetFillList.clear();
        for(int bpt : v_pt) {
          for(int bx : v_x) {
            for(int bz : v_z) {
              for(int bq : v_q) {
                for(int by : v_y) {
                  if(!CheckDiagonal(bpt,bx,bz)) {
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
  outfile->Close();
  cout << outfileN << " written." << endl;

  // call draw program
  /*
  TString cmd = "./draw.exe "+outfileN;
  system(cmd.Data());
  */

};

////////////////////////////////////////////////////
 
// return true, if `diagonalBinsOnly` mode is on, and this is an 
// off-diagonal bin
Bool_t CheckDiagonal(int cpt, int cx, int cz) {
  return diagonalBinsOnly &&
      ( cpt!=cx || cx!=cz );
};

// scan through array of bins `binArr`, checking each one; the vector `v` will
// contain the list of array indices for which the cut on `var` is satisfied
void CheckBins(TObjArray *binArr, std::vector<int> &v, Double_t var) {
  v.clear();
  for(int b=0; b<binArr->GetEntries(); b++) {
    if(((CutDef*)binArr->At(b))->Cut(var)) v.push_back(b);
  };
};
