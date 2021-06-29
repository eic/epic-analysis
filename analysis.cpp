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
  TObjArrayIter itParticle(tr->UseBranch("Particle"));
  TObjArrayIter itEFlowTrack(tr->UseBranch("EFlowTrack"));
  TObjArrayIter itEFlowPhoton(tr->UseBranch("EFlowPhoton"));
  TObjArrayIter itEFlowNeutralHadron(tr->UseBranch("EFlowNeutralHadron"));
  TObjArrayIter itPIDSystemsTrack(tr->UseBranch("PIDSystemsTrack"));
  
  
  // vars
  Double_t eleP,maxEleP;


  // define particle sets for histograms
  enum setEnum { pipTrack, nSets };
  Histos *histSet[nSets];
  histSet[pipTrack] = new Histos("pipTrack","#pi^{+} track");

  // define kinematics
  Kinematics *kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);

  // tree loop
  //ENT = 10000; // limiter
  Int_t s;
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

      // get hadronic final state variables
      kin->GetHadronicFinalState(itTrack, itEFlowTrack, itEFlowPhoton, itEFlowNeutralHadron, itPIDSystemsTrack, itParticle);
      
      // calculate DIS kinematics
      kin->CalculateDISbyElectron();



      // track loop
      itTrack.Reset();
      while(Track *trk = (Track*) itTrack()) {
        //cout << e << " " << trk->PID << endl;
        switch(trk->PID) {
          case 211: s=pipTrack; break;
          default: s=-1;
        };

        if(s>=0) {
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
	  
	  
	  
          // apply cuts and fill histograms
          if(kin->CutFull()) {
            // DIS kinematics
            histSet[s]->Hist("Q2vsX")->Fill(kin->x,kin->Q2);
            histSet[s]->Hist("W")->Fill(kin->W);
            histSet[s]->Hist("y")->Fill(kin->y);
            // hadron 4-momentum
            histSet[s]->Hist("p")->Fill(trk->P);
            histSet[s]->Hist("pTlab")->Fill(trk->PT);
            histSet[s]->Hist("eta")->Fill(trk->Eta);
            histSet[s]->Hist("phi")->Fill(trk->Phi);
            // hadron kinematics
            histSet[s]->Hist("z")->Fill(kin->z);
            histSet[s]->Hist("pT")->Fill(kin->pT);
            histSet[s]->Hist("qT")->Fill(kin->qT);
            histSet[s]->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2));
            histSet[s]->Hist("mX")->Fill(kin->mX);
            histSet[s]->Hist("phiH")->Fill(kin->phiH);
            histSet[s]->Hist("phiS")->Fill(kin->phiS);
          };
        };
      };
  };

  // write histograms
  outfile->cd();
  for(Int_t hs=0; hs<nSets; hs++) histSet[hs]->WriteHists();

};
