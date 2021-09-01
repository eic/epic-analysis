#include "AnalysisDD4hep.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

AnalysisDD4hep::AnalysisDD4hep(
  TString infileName_,
  Double_t eleBeamEn_,
  Double_t ionBeamEn_,
  Double_t crossingAngle_,
  TString outfilePrefix_
)
  : infileName(infileName_)
  , eleBeamEn(eleBeamEn_)
  , ionBeamEn(ionBeamEn_)
  , crossingAngle(crossingAngle_)
  , outfilePrefix(outfilePrefix_)  
{

  // set bin schemes
  AddBinScheme("pt","p_{T}");
  AddBinScheme("z","z");
  AddBinScheme("x","x");
  AddBinScheme("q","Q");
  AddBinScheme("y","y");
  AddBinScheme("pt_jet", "jet p_{T}");
  AddBinScheme("z_jet", "jet z");
  // final state bins (e.g., tracks or jets)
  AddBinScheme("finalState","finalState");
  AddFinalState("pipTrack","#pi^{+} tracks", 211);
  //AddFinalState("pimTrack","#pi^{-} tracks",-211);
  AddBinScheme("recMethod", "recMethod");
  AddRecMethod("Ele", "electron method");
  AddRecMethod("DA", "DA method");
  AddRecMethod("JB", "JB method");  
  // initialize diagonalizer settings
  // TODO: generalized diagonalizer
  diagonalPtXZ = false;
  diagonalXZQ = false;
  writeSimpleTree = false;
  maxEvents = 0;
};

// destructor
AnalysisDD4hep::~AnalysisDD4hep() {
};

void AnalysisDD4hep::process_event()
{
  cout << "-- running analysis of " << infileName << endl;

  // instantiate objects
  kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);
  kinTrue = new Kinematics(eleBeamEn, ionBeamEn, crossingAngle);
  ST = new SimpleTree("tree",kin);
  weight = new WeightsUniform();
  weightJet = new WeightsUniform();

  TChain *chain = new TChain("events");
  chain->Add(infileName);

  // FIXME: replace it with ExRootTreeReader::UseBranch()?
  TTreeReader tr(chain);

  // Truth
  TTreeReaderArray<Int_t>    mcparticles2_pdgID(tr, "mcparticles2.pdgID");
  TTreeReaderArray<Double_t> mcparticles2_psx(tr, "mcparticles2.psx");
  TTreeReaderArray<Double_t> mcparticles2_psy(tr, "mcparticles2.psy");
  TTreeReaderArray<Double_t> mcparticles2_psz(tr, "mcparticles2.psz");
  TTreeReaderArray<Int_t>    mcparticles2_status(tr, "mcparticles2.status");
  TTreeReaderArray<Int_t>    mcparticles2_genStatus(tr, "mcparticles2.genStatus");
  TTreeReaderArray<Double_t> mcparticles2_mass(tr, "mcparticles2.mass");

  // Reco
  TTreeReaderArray<Long64_t> ReconstructedParticles_pid(tr, "ReconstructedParticles.pid");
  TTreeReaderArray<Double_t> ReconstructedParticles_energy(tr, "ReconstructedParticles.energy");
  TTreeReaderArray<Double_t> ReconstructedParticles_p_x(tr, "ReconstructedParticles.p.x");
  TTreeReaderArray<Double_t> ReconstructedParticles_p_y(tr, "ReconstructedParticles.p.y");
  TTreeReaderArray<Double_t> ReconstructedParticles_p_z(tr, "ReconstructedParticles.p.z");
  TTreeReaderArray<Double_t> ReconstructedParticles_mass(tr, "ReconstructedParticles.mass");

  //HcalEndcap
  TTreeReaderArray<Float_t>  HcalHadronEndcapClusters_energy(tr, "HcalHadronEndcapClusters.energy");
  TTreeReaderArray<Float_t>  HcalHadronEndcapClusters_edep(tr, "HcalHadronEndcapClusters.edep");
  TTreeReaderArray<Double_t> HcalHadronEndcapClusters_x(tr, "HcalHadronEndcapClusters.position.x");
  TTreeReaderArray<Double_t> HcalHadronEndcapClusters_y(tr, "HcalHadronEndcapClusters.position.y");
  TTreeReaderArray<Double_t> HcalHadronEndcapClusters_z(tr, "HcalHadronEndcapClusters.position.z");
  TTreeReaderArray<Float_t>  HcalHadronEndcapClusters_theta(tr, "HcalHadronEndcapClusters.polar.theta");
  TTreeReaderArray<Float_t>  HcalHadronEndcapClusters_phi(tr, "HcalHadronEndcapClusters.polar.phi");

  //HcalBarrel
  TTreeReaderArray<Float_t>  HcalBarrelClusters_energy(tr, "HcalBarrelClusters.energy");
  TTreeReaderArray<Float_t>  HcalBarrelClusters_edep(tr, "HcalBarrelClusters.edep");
  TTreeReaderArray<Double_t> HcalBarrelClusters_x(tr, "HcalBarrelClusters.position.x");
  TTreeReaderArray<Double_t> HcalBarrelClusters_y(tr, "HcalBarrelClusters.position.y");
  TTreeReaderArray<Double_t> HcalBarrelClusters_z(tr, "HcalBarrelClusters.position.z");
  TTreeReaderArray<Float_t>  HcalBarrelClusters_theta(tr, "HcalBarrelClusters.polar.theta");
  TTreeReaderArray<Float_t>  HcalBarrelClusters_phi(tr, "HcalBarrelClusters.polar.phi");

  //Ecal
  TTreeReaderArray<Float_t>  EcalEndcapPClusters_energy(tr, "EcalEndcapPClusters.energy");
  TTreeReaderArray<Float_t>  EcalEndcapPClusters_edep(tr, "EcalEndcapPClusters.edep");
  TTreeReaderArray<Double_t> EcalEndcapPClusters_x(tr, "EcalEndcapPClusters.position.x");
  TTreeReaderArray<Double_t> EcalEndcapPClusters_y(tr, "EcalEndcapPClusters.position.y");
  TTreeReaderArray<Double_t> EcalEndcapPClusters_z(tr, "EcalEndcapPClusters.position.z");
  TTreeReaderArray<Float_t>  EcalEndcapPClusters_theta(tr, "EcalEndcapPClusters.polar.theta");
  TTreeReaderArray<Float_t>  EcalEndcapPClusters_phi(tr, "EcalEndcapPClusters.polar.phi");

  TTreeReaderArray<Float_t>  EcalEndcapNClusters_energy(tr, "EcalEndcapNClusters.energy");
  TTreeReaderArray<Float_t>  EcalEndcapNClusters_edep(tr, "EcalEndcapNClusters.edep");
  TTreeReaderArray<Double_t> EcalEndcapNClusters_x(tr, "EcalEndcapNClusters.position.x");
  TTreeReaderArray<Double_t> EcalEndcapNClusters_y(tr, "EcalEndcapNClusters.position.y");
  TTreeReaderArray<Double_t> EcalEndcapNClusters_z(tr, "EcalEndcapNClusters.position.z");
  TTreeReaderArray<Float_t>  EcalEndcapNClusters_theta(tr, "EcalEndcapNClusters.polar.theta");
  TTreeReaderArray<Float_t>  EcalEndcapNClusters_phi(tr, "EcalEndcapNClusters.polar.phi");

  TTreeReaderArray<Float_t>  EcalBarrelClusters_energy(tr, "EcalBarrelClusters.energy");
  TTreeReaderArray<Float_t>  EcalBarrelClusters_edep(tr, "EcalBarrelClusters.edep");
  TTreeReaderArray<Double_t> EcalBarrelClusters_x(tr, "EcalBarrelClusters.position.x");
  TTreeReaderArray<Double_t> EcalBarrelClusters_y(tr, "EcalBarrelClusters.position.y");
  TTreeReaderArray<Double_t> EcalBarrelClusters_z(tr, "EcalBarrelClusters.position.z");
  TTreeReaderArray<Float_t>  EcalBarrelClusters_theta(tr, "EcalBarrelClusters.polar.theta");
  TTreeReaderArray<Float_t>  EcalBarrelClusters_phi(tr, "EcalBarrelClusters.polar.phi");

  TTreeReader::EEntryStatus entrystats = tr.SetEntry(0);

  // open output file
  // FIXME: replace name variable
  TFile *outfile = new TFile(outfilePrefix,"RECREATE");

  // number of bins
  const Int_t NptBins = BinScheme("pt")->GetNumBins();
  const Int_t NxBins = BinScheme("x")->GetNumBins();
  const Int_t NzBins = BinScheme("z")->GetNumBins();
  const Int_t NqBins = BinScheme("q")->GetNumBins();
  const Int_t NyBins = BinScheme("y")->GetNumBins();
  const Int_t NfinalStateBins = BinScheme("finalState")->GetNumBins();
  const Int_t NptjetBins = BinScheme("pt_jet")->GetNumBins();
  const Int_t NzjetBins = BinScheme("z_jet")->GetNumBins();
  const Int_t NrecMethodBins = BinScheme("recMethod")->GetNumBins();

  // sets of histogram sets
  // - `histSet` is a data structure for storing and organizing pointers to
  //   sets of histograms (`Histos` objects)
  // - `histSet*List` are used as temporary lists of relevant `Histos` pointers
  // - TODO: if we add one more dimension, 7D array will probably break; need
  //         better data structure
  Histos *histSet[NptBins][NxBins][NzBins][NqBins][NyBins][NfinalStateBins];
  Histos *histSetJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins];
  Histos *histSetBreitJets[NptjetBins][NzjetBins][NxBins][NqBins][NyBins][NrecMethodBins];

  std::vector<Histos*> histSetList;
  std::vector<Histos*> histSetListJets;
  std::vector<Histos*> histSetListBreitJets;

  std::vector<Histos*> histSetFillList;
  std::vector<int> v_pt, v_x, v_z, v_q, v_y;
  // instantiate Histos sets, and populate 
  TString histosN,histosT;

  cout << "Define track histograms..." << endl;
  for(int bpt=0; bpt<NptBins; bpt++) { // - loop over pT bins
    for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
      for(int bz=0; bz<NzBins; bz++) { // - loop over z bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q bins
          if(CheckDiagonal(bpt,bx,bz,bq)) continue; // diagonalizer
          for(int by=0; by<NyBins; by++) { // - loop over y bins
            for(int bfs=0; bfs<NfinalStateBins; bfs++) { // - loop over final states

              // set Histos name and title
              histosN = this->GetHistosName (bpt,bx,bz,bq,by,bfs);
              histosT = this->GetHistosTitle(bpt,bx,bz,bq,by,bfs);

              // define set of histograms for this bin
              histSet[bpt][bx][bz][bq][by][bfs] = new Histos(histosN,histosT);
              HS = histSet[bpt][bx][bz][bq][by][bfs]; // shorthand pointer

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
              // -- DIS kinematics resolution
              HS->DefineHist1D("xRes","x - x_{true}","", NBINS, -1, 1);
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

              // store cut definitions with histogram sets
              HS->AddCutDef(BinScheme("pt")->Cut(bpt));
              HS->AddCutDef(BinScheme("x")->Cut(bx));
              HS->AddCutDef(BinScheme("z")->Cut(bz));
              HS->AddCutDef(BinScheme("q")->Cut(bq));
              HS->AddCutDef(BinScheme("y")->Cut(by));
              HS->AddCutDef(BinScheme("finalState")->Cut(bfs));

              // add histogram set full list
              histSetList.push_back(histSet[bpt][bx][bz][bq][by][bfs]);
            };
          };
        };
      };
    };
  };

  cout << "Define jet histograms..." << endl;
  for(int bpt=0; bpt<NptjetBins; bpt++) { // - loop over jet pT bins
    for(int bz=0; bz<NzjetBins; bz++){
      for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q2 bins
          for(int by=0; by<NyBins; by++) { // - loop over y bins

            histosN = this->GetHistosNameJets(bpt, bz, bx, bq, by);
            histosT = this->GetHistosTitleJets(bpt, bz, bx, bq, by);

            histSetJets[bpt][bz][bx][bq][by] = new Histos(histosN,histosT);
            HS = histSetJets[bpt][bz][bx][bq][by]; // shorthand pointer

            // jet kinematics plots
            HS->DefineHist1D("pT_jet","p_{T}","GeV", NBINS, 1e-2, 50);
            HS->DefineHist1D("mT_jet","m_{T}","GeV", NBINS, 1e-2, 20);
            HS->DefineHist1D("z_jet","z","GeV", NBINS,0, 1);
            HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
            HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0);
            HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
            HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);
            // store cut definitions with histogram sets, then add histogram sets full list
            HS->AddCutDef(BinScheme("pt_jet")->Cut(bpt));
            HS->AddCutDef(BinScheme("z_jet")->Cut(bz));
            HS->AddCutDef(BinScheme("x")->Cut(bx));
            HS->AddCutDef(BinScheme("q2")->Cut(bq));
            HS->AddCutDef(BinScheme("y")->Cut(by));
            histSetListJets.push_back(histSetJets[bpt][bz][bx][bq][by]);
          };	
        };
      };
    };
  };
  #if INCCENTAURO == 1
  for(int bpt=0; bpt<NptjetBins; bpt++) { // - loop over jet pT bins
    for(int bz=0; bz<NzjetBins; bz++){
      for(int bx=0; bx<NxBins; bx++) { // - loop over x bins
        for(int bq=0; bq<NqBins; bq++) { // - loop over q2 bins
          for(int by=0; by<NyBins; by++) { // - loop over y bins
            for(int brec=0; brec<NrecMethodBins; brec++){
              histosN = this->GetHistosNameBreitJets(bpt, bz, bx, bq, by, brec);
              histosT = this->GetHistosTitleBreitJets(bpt, bz, bx, bq, by, brec);
              histSetBreitJets[bpt][bz][bx][bq][by][brec] = new Histos(histosN,histosT);
              HS = histSetBreitJets[bpt][bz][bx][bq][by][brec]; // shorthand pointer

              // jet kinematics plots
              HS->DefineHist1D("pT_jet","p_{T}","GeV", NBINS, 1e-2, 20);
              HS->DefineHist1D("mT_jet","m_{T}","GeV", NBINS, 1e-2, 20);
              HS->DefineHist1D("z_jet","z","GeV", NBINS,0, 1);
              HS->DefineHist1D("eta_jet","#eta_{lab}","GeV", NBINS,-5,5);
              HS->DefineHist1D("qT_jet","qT", "GeV", NBINS, 0, 10.0);
              HS->DefineHist1D("jperp","j_{#perp}","GeV", NBINS, 0, 3.0);
              HS->DefineHist1D("qTQ", "q_{T}/Q, jets", "GeV", NBINS, 0, 3.0);

              // store cut definitions with histogram sets, then add histogram sets full list
              HS->AddCutDef(BinScheme("pt_jet")->Cut(bpt));
              HS->AddCutDef(BinScheme("z_jet")->Cut(bz));
              HS->AddCutDef(BinScheme("x")->Cut(bx));
              HS->AddCutDef(BinScheme("q2")->Cut(bq));
              HS->AddCutDef(BinScheme("y")->Cut(by));
              HS->AddCutDef(BinScheme("recMethod")->Cut(brec));
              histSetListBreitJets.push_back(histSetBreitJets[bpt][bz][bx][bq][by][brec]);
            };
          };
        };
      };
    };
  };
  #endif

  
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
  Long64_t numGen = tr.GetEntries();
  TString sep = "--------------------------------------------";
  cout << sep << endl;
  cout << "assumed total cross section: " << xsecTot << " nb" << endl;
  cout << "number of generated events:  " << numGen << endl;

  // initialize total weights
  Double_t wTotal = 0.;
  Double_t wJetTotal = 0.;

  // vars
  Double_t eleP,maxEleP;
  int pid,bFinalState;
  Double_t elePtrue, maxElePtrue;

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  while(tr.Next())
    {
      
      double maxP = 0;
      for(int imc=0; imc>mcparticles2_pdgID.GetSize(); imc++)
	{
	  // genStatus 4: beam particle 1: final state 
	  if(mcparticles2_pdgID[imc] == 11 && mcparticles2_genStatus[imc] == 1)
	    {
	      double p_ = sqrt(pow(mcparticles2_psx[imc],2) + pow(mcparticles2_psy[imc],2) + pow(mcparticles2_psz[imc],2));
	      double mass_ = mcparticles2_mass[imc]; // in GeV
	      if(p_ > maxP)
		{
		  maxP = p_;
		  kinTrue->vecElectron.SetPxPyPzE(mcparticles2_psx[imc],
						  mcparticles2_psy[imc],
						  mcparticles2_psz[imc],
						  sqrt(p_*p_ + mass_*mass_));
		}
	    }// if electron
	}//mcparticles loop
      
      // Loop over calorimeters
      // fill cluster container
      vector<Clusters*> v_ecal_clusters;
      vector<Clusters*> v_hcal_clusters;

      // HcalBarrel
      for(int icl=0; icl<HcalBarrelClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = HcalBarrelClusters_energy[icl]; // use edep?
	  clus->x = HcalBarrelClusters_x[icl]; // use edep?
	  clus->y = HcalBarrelClusters_y[icl]; // use edep?
	  clus->z = HcalBarrelClusters_z[icl]; // use edep?
	  clus->theta = HcalBarrelClusters_theta[icl]; // use edep?
	  clus->phi = HcalBarrelClusters_phi[icl]; // use edep?
	  v_hcal_clusters.push_back(clus);
	}
      // HcalEndcap
      for(int icl=0; icl<HcalHadronEndcapClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = HcalHadronEndcapClusters_energy[icl]; // use edep?
	  clus->x = HcalHadronEndcapClusters_x[icl]; // use edep?
	  clus->y = HcalHadronEndcapClusters_y[icl]; // use edep?
	  clus->z = HcalHadronEndcapClusters_z[icl]; // use edep?
	  clus->theta = HcalHadronEndcapClusters_theta[icl]; // use edep?
	  clus->phi = HcalHadronEndcapClusters_phi[icl]; // use edep?
	  v_hcal_clusters.push_back(clus);
	}

      // Barrel
      for(int icl=0; icl<EcalBarrelClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = EcalBarrelClusters_energy[icl]; // use edep?
	  clus->x = EcalBarrelClusters_x[icl]; // use edep?
	  clus->y = EcalBarrelClusters_y[icl]; // use edep?
	  clus->z = EcalBarrelClusters_z[icl]; // use edep?
	  clus->theta = EcalBarrelClusters_theta[icl]; // use edep?
	  clus->phi = EcalBarrelClusters_phi[icl]; // use edep?
	  v_ecal_clusters.push_back(clus);
	}

      // EndcapP
      for(int icl=0; icl<EcalEndcapPClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = EcalEndcapPClusters_energy[icl]; // use edep?
	  clus->x = EcalEndcapPClusters_x[icl]; // use edep?
	  clus->y = EcalEndcapPClusters_y[icl]; // use edep?
	  clus->z = EcalEndcapPClusters_z[icl]; // use edep?
	  clus->theta = EcalEndcapPClusters_theta[icl]; // use edep?
	  clus->phi = EcalEndcapPClusters_phi[icl]; // use edep?
	  v_ecal_clusters.push_back(clus);
	}

      // EndcapN
      for(int icl=0; icl<EcalEndcapNClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = EcalEndcapNClusters_energy[icl]; // use edep?
	  clus->x = EcalEndcapNClusters_x[icl]; // use edep?
	  clus->y = EcalEndcapNClusters_y[icl]; // use edep?
	  clus->z = EcalEndcapNClusters_z[icl]; // use edep?
	  clus->theta = EcalEndcapNClusters_theta[icl]; // use edep?
	  clus->phi = EcalEndcapNClusters_phi[icl]; // use edep?
	  v_ecal_clusters.push_back(clus);
	}

      // find scattered electron
      // FIXME: hard-coded e_threshold
      int electron_index = find_electron(v_ecal_clusters, v_hcal_clusters, 3.0);
      if(electron_index < 0)
	{
	  cout << "No scattered electron found" << endl;
	  continue;
	}

      double electron_E     = v_ecal_clusters[electron_index]->E/1000.0; 
      double electron_theta = v_ecal_clusters[electron_index]->theta;
      double electron_phi = v_ecal_clusters[electron_index]->phi;
      // FIXME: use track information?
      double electron_pt    = electron_E * sin(electron_theta);
      double electron_eta = -log(tan(electron_theta/2.0));

      kin->vecElectron.SetPtEtaPhiM(
				    electron_pt,
				    electron_eta,
				    electron_phi,
				    Kinematics::ElectronMass()
				    );

      // calculate DIS kinematics
      kin->CalculateDISbyElectron(); // reconstructed
      kinTrue->CalculateDISbyElectron(); // generated (truth)

      TLorentzVector v_had;
      double hpx=0; 
      double hpy=0;
      for(int itrk=0; itrk<ReconstructedParticles_pid.GetSize(); itrk++)
	{
	  // FIXME: I think pid is using the true information
	  // Add PID smearing here 
	  int pid = abs(ReconstructedParticles_pid[itrk]);

	  TLorentzVector v_temp;
	  double reco_E = ReconstructedParticles_energy[itrk];
	  double reco_px = ReconstructedParticles_p_x[itrk];
	  double reco_py = ReconstructedParticles_p_y[itrk];
	  double reco_pz = ReconstructedParticles_p_z[itrk];
	  double reco_mass = ReconstructedParticles_mass[itrk];
	  v_temp.SetPxPyPzE(reco_px, reco_py, reco_pz, reco_E);
	  kin->vecHadron.SetPxPyPzE(
				    reco_px,
				    reco_py,
				    reco_pz,
				    reco_E);
				      
	  kin->CalculateHadronKinematics();
	  //FIXME: add true information for all hadrons

	  Double_t w = weight->GetWeight(*kin);
	  wTotal += w;

	  hpx += reco_px;
	  hpy += reco_py;
	  v_had += v_temp;

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
	    CheckBins( BinScheme("q2"), v_q,  kin->Q2 );
	    CheckBins( BinScheme("y"),  v_y,  kin->y );

	    // build list of histogram sets to fill
	    histSetFillList.clear();
	    for(int bpt : v_pt) {
	      for(int bx : v_x) {
		for(int bz : v_z) {
		  for(int bq : v_q) {
		    for(int by : v_y) {
		      if(!CheckDiagonal(bpt,bx,bz,bq)) {
			histSetFillList.push_back(histSet[bpt][bx][bz][bq][by][bFinalState]);
		      };
		    };
		  };
		};
	      };
	    };

	    // loop through list of histogram sets, and fill them
	    for(Histos *H : histSetFillList) {
	      // DIS kinematics
	      dynamic_cast<TH2*>(H->Hist("Q2vsX"))->Fill(kin->x,kin->Q2,w);
	      H->Hist("Q")->Fill(TMath::Sqrt(kin->Q2),w);
	      H->Hist("x")->Fill(kin->x,w);
	      H->Hist("W")->Fill(kin->W,w);
	      H->Hist("y")->Fill(kin->y,w);
	      // hadron 4-momentum
	      H->Hist("pLab")->Fill(kin->pLab,w);
	      H->Hist("pTlab")->Fill(kin->pTlab,w);
	      H->Hist("etaLab")->Fill(kin->etaLab,w);
	      H->Hist("phiLab")->Fill(kin->phiLab,w);
	      // hadron kinematics
	      H->Hist("z")->Fill(kin->z,w);
	      H->Hist("pT")->Fill(kin->pT,w);
	      H->Hist("qT")->Fill(kin->qT,w);
	      H->Hist("qTq")->Fill(kin->qT/TMath::Sqrt(kin->Q2),w);
	      H->Hist("mX")->Fill(kin->mX,w);
	      H->Hist("phiH")->Fill(kin->phiH,w);
	      H->Hist("phiS")->Fill(kin->phiS,w);
	      // cross sections (divide by lumi after all events processed)
	      H->Hist("Q_xsec")->Fill(TMath::Sqrt(kin->Q2),w);
	      // DIS kinematics resolution
	      //	      H->Hist("xRes")->Fill(kin->x - kinTrue->x,w);
	    };

	    // fill simple tree (not binned)
	    // TODO: consider adding a `finalState` cut
	    if( writeSimpleTree && histSetFillList.size()>0 ) ST->FillTree(w);
	
	  }//if cut
	}//trk loop

      //Hadronic reconstruction method
      kin->sigmah = (kin->vecElectron.E() - kin->vecElectron.Pz());
      kin->Pxh = hpx - kin->vecElectron.Px();
      kin->Pyh = hpy - kin->vecElectron.Py();

    }// tree reader loop

}//execute

double AnalysisDD4hep::isolation(double cone_theta, double cone_phi, std::vector<Clusters*> cluster_container, double E_threshold)
{
  double cone_eta = -log(tan(cone_theta/2.0));
  double cone_iso = 0.0;

  // loop over clusters
  for(auto icl : cluster_container)
    {
      auto clusters = *icl;
      double clus_e     = clusters.E;
      double clus_phi   = clusters.phi;
      double clus_theta = clusters.theta;
      double clus_eta   = -log(tan(clus_theta/2.0));
      double clus_pt    = clus_e*sin(clus_theta);
      
      double dphi = clus_phi - cone_phi;
      //      dphi = (dphi + TMath::Pi()) % (2 * TMath::Pi()) - TMath::Pi();
      double dr = sqrt( pow(dphi,2) + pow((clus_eta - cone_eta), 2) );

      // get E_cone
      // FIXME: using hard-coded cone radius
      if(dr < 1.0)
	{
	  cone_iso += clus_e;
	}
    }// for loop

  return cone_iso;
}

int AnalysisDD4hep::find_electron(std::vector<Clusters*> ecal_cluster, std::vector<Clusters*> hcal_cluster, double e_threshold)
{

  // init 
  double ptmax = 0;
  double electron_iso = 999;
  int index_max = -999;
  
  int index = 0;
  for(auto icl : ecal_cluster)
    {
      index++;
      auto icluster = *icl;
      double clus_E = icluster.E/1000.0;
      if(clus_E < e_threshold)
	continue;
      if(icluster.theta * 180./TMath::Pi() < 2)
	continue;

      double clus_theta = icluster.theta;
      double clus_phi   = icluster.phi;
      double clus_eta   = -log(tan(clus_theta/2.0));
      double clus_pt    = clus_E * sin(clus_theta);
      
      // FIXME hard-coded threshold 
      double clus_ecal_iso = isolation(clus_theta, clus_phi, ecal_cluster, 0.1);
      double clus_hcal_iso = isolation(clus_theta, clus_phi, hcal_cluster, 0.1);
      double clus_iso = clus_ecal_iso + clus_hcal_iso - clus_E; // subtract the electron energy

      // FIXME: using hard-coded isolation criteria; 10%
      if(clus_iso > clus_E*0.1)
	continue;

      if(clus_pt > ptmax)
	{
	  ptmax = clus_pt;
	  electron_iso = clus_iso;
	  index_max = index;
	}
      }// loop over ecal clusters
  
  return index_max-1;
}
