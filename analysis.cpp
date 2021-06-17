// root
#include "TChain.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TFile.h"

// delphes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

// largex-eic
#include "Histos.h"
#include "Kinematics.h"


//=========================================================================

int main(int argc, char **argv) {

  // ARGUMENTS ////////////////////////////////////////////////
  TString infile="datarec/test_crossDivNrgCrab_25mRad_5x41_v1.root";
  Double_t eleBeamEn = 5; // GeV
  Double_t ionBeamEn = 41; // GeV
  Double_t crossingAngle = 25; // mrad
  if(argc>1) infile = TString(argv[1]);
  if(argc>2) eleBeamEn = (Double_t)strtof(argv[2],NULL);
  if(argc>3) ionBeamEn = (Double_t)strtof(argv[3],NULL);
  /////////////////////////////////////////////////////////////


  // read delphes tree
  TChain *chain = new TChain("Delphes");
  chain->Add(infile);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  Long64_t ENT = tr->GetEntries();

  // define output file
  TFile *outfile = new TFile("histos.root","RECREATE");

  // branch iterators
  TObjArrayIter itTrack(tr->UseBranch("Track"));
  TObjArrayIter itElectron(tr->UseBranch("Electron"));

  // vars
  Double_t eleP,maxEleP;


  // define particle sets for histograms
  enum setEnum { pipTrack, nSets };
  Histos *histSet[nSets];
  histSet[pipTrack] = new Histos("pipTrack","#pi^{+} track");

  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);

  // tree loop
  ENT = 1000; // limiter
  Int_t s;
  for(Long64_t e=0; e<ENT; e++) {
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
      kin->DISbyElectron();



      // track loop
      itTrack.Reset();
      while(Track *trk = (Track*) itTrack()) {
        cout << e << " " << trk->PID << endl;
        switch(trk->PID) {
          case 211: s=pipTrack; break;
          default: s=-1;
        };

        if(s>=0) {
          histSet[s]->Hist("p")->Fill(trk->P);
          histSet[s]->Hist("pT")->Fill(trk->PT);
          histSet[s]->Hist("eta")->Fill(trk->Eta);
          histSet[s]->Hist("phi")->Fill(trk->Phi);
        };
      };
  };

  // write histograms
  outfile->cd();
  for(Int_t hs=0; hs<nSets; hs++) histSet[hs]->WriteHists();

};
