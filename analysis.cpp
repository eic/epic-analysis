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
  // - y-minimum cuts
  const int NY=4;
  Float_t ycut[NY] = {
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
  int s_pid,pid;

  // build array of histogram sets
  Histos *histSet[NY][NPart];
  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetFillList;
  TString keyN,keyT;
  // - loop over y cuts
  for(int hy=0; hy<NY; hy++) {
    keyN = Form("_ycut%d",hy);
    keyT = Form(", y>%.2f",ycut[hy]);
    // loop over particles
    histSet[hy][pPip] = new Histos("pipTrack"+keyN,"#pi^{+} tracks"+keyT);
    histSet[hy][pPim] = new Histos("pimTrack"+keyN,"#pi^{-} tracks"+keyT);
    // add to full list
    for(int hp=0; hp<NPart; hp++) histSetList.push_back(histSet[hy][hp]);
  };

  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);

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

        // chosen PID value
        pid = trk->PID;

        // decide which histogram sets to fill
        // - check PID, to see if it's a particle we're interested in for
        //   histograms; if not, proceed to next
        auto kv = PIDtoEnum.find(pid);
        if(kv!=PIDtoEnum.end()) s_pid = kv->second;
        else continue;
        // - check y cut; `v_y` will be the list of `s_y`s for which
        //   the y cut is satisfied
        v_y.clear();
        for(int s_y=0; s_y<NY; s_y++) {
          if( kin->y > ycut[s_y] ) v_y.push_back(s_y);
        };

        // build list of histogram sets to fill, and proceed only 
        // if there are some sets
        histSetFillList.clear();
        for(int s_y : v_y) {
          histSetFillList.push_back(histSet[s_y][s_pid]);
        };
        if(histSetFillList.size()>0) {

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
            };
          };
        };
      };
  };

  // write histograms
  outfile->cd();
  for(Histos *H : histSetList) H->WriteHists();
  for(Histos *H : histSetList) H->Write();
  outfile->Close();
  cout << outfileN << " written." << endl;

  // call draw program
  TString cmd = "./draw.exe "+outfileN;
  system(cmd.Data());

};
