#include <stdlib.h>

// root
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TRegexp.h"

// delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

// largex-eic
#include "Histos.h"
#include "Kinematics.h"

// subroutines
void CenterDelta(Double_t center, Double_t delta, Double_t &cutmin, Double_t &cutmax);
Bool_t CheckDiagonal(int cpt, int cx, int cz);

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
  
  // ----------------------------------------------------------
  // BINNING
  // ----------------------------------------------------------

  // - pT bins
  const int NbinsPT = 4;
  Double_t ptmin[NbinsPT], ptmax[NbinsPT];
  ptmin[0]=0.; ptmax[0]=1000.; // element 0 is always full range
  CenterDelta( 0.15, 0.05, ptmin[1], ptmax[1]); // slide 14
  CenterDelta( 0.55, 0.05, ptmin[2], ptmax[2]); // slide 13
  CenterDelta( 0.50, 0.05, ptmin[3], ptmax[3]); // slide 11

  // - x bins
  const int NbinsX = 4;
  Double_t xmin[NbinsX], xmax[NbinsX];
  xmin[0]=0.; xmax[0]=1.; // element 0 is always full range
  CenterDelta( 0.1, 0.05, xmin[1], xmax[1]);
  CenterDelta( 0.6, 0.05, xmin[2], xmax[2]);
  CenterDelta( 0.3, 0.05, xmin[3], xmax[3]);

  // - z bins
  const int NbinsZ = 4;
  Double_t zmin[NbinsZ], zmax[NbinsZ];
  zmin[0]=0.; zmax[0]=1.; // element 0 is always full range
  CenterDelta( 0.7, 0.05, zmin[1], zmax[1] );
  CenterDelta( 0.5, 0.05, zmin[2], zmax[2] );
  CenterDelta( 0.7, 0.05, zmin[3], zmax[3] );

  // if this is true, only take 'diagonal' elements of the multi
  // dimensional array of possible pT,x,z bins; this is useful
  // if you want to check specific bins
  diagonalBinsOnly = true;


  // - y-minimum cuts
  const int NY=1;
  Double_t ymin[NY] = {
    0.00
    //0.03,
    //0.05,
    //0.10
  };

  // - particle species
  std::map<int,int> PIDtoEnum;
  enum partEnum{
    pPip,
    //pPim,
    NPart
  };
  PIDtoEnum.insert(std::pair<int,int>(211,pPip));
  //PIDtoEnum.insert(std::pair<int,int>(-211,pPim));

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



  // sets of histogram sets
  // - `histSet` is a data structure for storing and organizing pointers to
  //   sets of histograms (`Histos` objects)
  // - `histSet*List` are used as temporary lists of relevant `Histos` pointers
  Histos *histSet[NbinsPT][NbinsX][NbinsZ][NY][NPart];
  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetFillList;
  std::vector<int> v_pt, v_x, v_z, v_y;
  // instantiate Histos sets, and populate 
  TString keyN,keyT;
  cout << "Define histograms..." << endl;
  for(int bpt=0; bpt<NbinsPT; bpt++) { // - loop over pT bins
    for(int bx=0; bx<NbinsX; bx++) { // - loop over x bins
      for(int bz=0; bz<NbinsZ; bz++) { // - loop over z bins
        if(CheckDiagonal(bpt,bx,bz)) continue;
        for(int by=0; by<NY; by++) { // - loop over y cuts
          // set plot name
          keyN  = Form("_ptbin%d",bpt);
          keyN += Form("_xbin%d",bx);
          keyN += Form("_zbin%d",bz);
          keyN += Form("_ycut%d",by);
          // set plot title
          keyT  = bpt>0 ? Form(", %.2f<p_{T}<%.2f",ptmin[bpt],ptmax[bpt])
                        : ", full p_{T}";
          keyT += bx>0 ? Form(", %.2f<x<%.2f",xmin[bx],xmax[bx])
                       : ", full x";
          keyT += bz>0 ? Form(", %.2f<z<%.2f",zmin[bz],zmax[bz])
                       : ", full z";
          keyT += Form(", y>%.2f",ymin[by]);
          // loop over particles
          histSet[bpt][bx][bz][by][pPip] = new Histos("pipTrack"+keyN,"#pi^{+} tracks"+keyT);
          //histSet[bpt][bx][bz][by][pPim] = new Histos("pimTrack"+keyN,"#pi^{-} tracks"+keyT);
          // add to full list
          for(int bp=0; bp<NPart; bp++) histSetList.push_back(histSet[bpt][bx][bz][by][bp]);
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
        v_pt.clear();
        for(int bpt=0; bpt<NbinsPT; bpt++) {
          if( kin->pT>ptmin[bpt] && kin->pT<ptmax[bpt] ) v_pt.push_back(bpt);
        };
        // - check x bin
        v_x.clear();
        for(int bx=0; bx<NbinsX; bx++) {
          if( kin->x>xmin[bx] && kin->x<xmax[bx] ) v_x.push_back(bx);
        };
        // - check z bin
        v_z.clear();
        for(int bz=0; bz<NbinsZ; bz++) {
          if( kin->z>zmin[bz] && kin->z<zmax[bz] ) v_z.push_back(bz);
        };
        // - check y cut
        v_y.clear();
        for(int by=0; by<NY; by++) {
          if( kin->y > ymin[by] ) v_y.push_back(by);
        };

        // build list of histogram sets to fill
        histSetFillList.clear();
        for(int bpt : v_pt) {
          for(int bx : v_x) {
            for(int bz : v_z) {
              for(int by : v_y) {
                if(!CheckDiagonal(bpt,bx,bz)) {
                  histSetFillList.push_back(histSet[bpt][bx][bz][by][bpart]);
                };
              };
            };
          };
        };

        // loop through list of histogram sets, and fill them
        for(Histos *H : histSetFillList) {
          // DIS kinematics
          H->Hist("Q2vsX")->Fill(kin->x,kin->Q2);
          H->Hist("W")->Fill(kin->W);
          H->Hist("y")->Fill(kin->y);
          // hadron 4-momentum
          H->Hist("p")->Fill(trk->P);
          H->Hist("pTlab")->Fill(trk->PT);
          H->Hist("eta")->Fill(trk->Eta);
          H->Hist("phi")->Fill(trk->Phi);
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
 
// define cut using `center` +/- `delta`; set cut minimum and maximum
void CenterDelta(Double_t center, Double_t delta, Double_t &cutmin, Double_t &cutmax) {
  cutmin = center - delta;
  cutmax = center + delta;
};

// return true, if `diagonalBinsOnly` mode is on, and this is an 
// off-diagonal bin
Bool_t CheckDiagonal(int cpt, int cx, int cz) {
  return diagonalBinsOnly &&
      ( cpt!=cx || cx!=cz );
};
