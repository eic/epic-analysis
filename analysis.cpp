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


//=========================================================================

int main(int argc, char **argv) {

  // ARGUMENTS ////////////////////////////////////////////////
  TString infile="datarec/test_crossDivNrgCrab_25mRad_5x41_v1.root";
  if(argc>1) infile = TString(argv[1]);
  /////////////////////////////////////////////////////////////


  // read delphes tree
  TChain *chain = new TChain("Delphes");
  chain->Add(infile);
  ExRootTreeReader *tr = new ExRootTreeReader(chain);
  Long64_t ENT = tr->GetEntries();

  // define output file
  TFile *outfile = new TFile("histos.root","RECREATE");

  // branches
  TClonesArray *brTrack = tr->UseBranch("Track");
  TObjArrayIter iterTrack(brTrack);

  // objects
  Track *trk;

  // define particle sets for histograms
  enum setEnum { pipTrack, nSets };
  Histos *histSet[nSets];
  histSet[pipTrack] = new Histos("pipTrack","#pi^{+} track");

  // tree loop
  ENT = 1000; // limiter
  Int_t s;
  for(Long64_t e=0; e<ENT; e++) {
    tr->ReadEntry(e);

      iterTrack.Reset();
      while((trk = (Track*) iterTrack())) {
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
