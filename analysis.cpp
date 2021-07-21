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
void DrawRatios(Histos *numerSet, Histos *denomSet);

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

  // vars
  Double_t eleP,maxEleP;


  // define histogram sets
  // - pT bins
  const int NbinsPT = 2;
  Double_t ptmin[NbinsPT];
  Double_t ptmax[NbinsPT];
  ptmin[0]=0.;         ptmax[0]=1000.;    // element 0 is always full range
  ptmin[1]=0.5-0.01;   ptmax[1]=0.5+0.01; // test bin
  std::vector<int> v_pt;
  // - x bins
  const int NbinsX = 2;
  Double_t xmin[NbinsX];
  Double_t xmax[NbinsX];
  xmin[0]=0.;          xmax[0]=1.;       // element 0 is always full range
  xmin[1]=0.3-0.01;    xmax[1]=0.3+0.01; // test bin
  std::vector<int> v_x;
  // - z bins
  const int NbinsZ = 2;
  Double_t zmin[NbinsZ];
  Double_t zmax[NbinsZ];
  zmin[0]=0.;          zmax[0]=1.;       // element 0 is always full range
  zmin[1]=0.3-0.01;    zmax[1]=0.3+0.01; // test bin
  std::vector<int> v_z;
  // - y-minimum cuts
  const int NY=4;
  Double_t ymin[NY] = {
    0.00,
    0.03,
    0.05,
    0.10
  };
  std::vector<int> v_y;
  // - particle species
  enum partEnum{
    pPip,
    pPim,
    NPart
  };
  std::map<int,int> PIDtoEnum;
  PIDtoEnum.insert(std::pair<int,int>(211,pPip));
  PIDtoEnum.insert(std::pair<int,int>(-211,pPim));
  int pid,bpart;


  // build array of histogram sets
  Histos *histSet[NbinsPT][NbinsX][NbinsZ][NY][NPart];
  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetFillList;
  TString keyN,keyT;
  for(int bpt=0; bpt<NbinsPT; bpt++) { // - loop over pT bins
    for(int bx=0; bx<NbinsX; bx++) { // - loop over x bins
      for(int bz=0; bz<NbinsZ; bz++) { // - loop over z bins
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
          histSet[bpt][bx][bz][by][pPim] = new Histos("pimTrack"+keyN,"#pi^{-} tracks"+keyT);
          // add to full list
          for(int bp=0; bp<NPart; bp++) histSetList.push_back(histSet[bpt][bx][bz][by][bp]);
        };
      };
    };
  };

  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);

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
  Double_t lumi = tr->GetEntries() / xsecTot; // [nb^-1]

  // tree loop
  //ENT = 1000; // limiter
  for(Long64_t e=0; e<ENT; e++) {
    if(e>0&&e%1000==0) cout << (Double_t)e/ENT*100 << "%" << endl;
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
                histSetFillList.push_back(histSet[bpt][bx][bz][by][bpart]);
              };
            };
          };
        };

        // loop through list of histogram sets, and fill them
        if(histSetFillList.size()>0) {
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
  };


  // write histograms
  outfile->cd();
  for(Histos *H : histSetList) H->WriteHists(outfile);
  for(Histos *H : histSetList) H->Write();
  outfile->Close();
  cout << outfileN << " written." << endl;

  // call draw program
  TString cmd = "./draw.exe "+outfileN;
  system(cmd.Data());

};
