#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include "AnalysisDD4hep.h"

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
) : Analysis(
  infileName_,
  eleBeamEn_,
  ionBeamEn_,
  crossingAngle_,
  outfilePrefix_
) {
  // dd4hep-specific settings defaults
  // initialize scatt. electron cuts
  fEThreshold = eleBeamEn_*0.1; // min energy cut
  fIsoR = 1.0;                  // Isolation cone R
  fIsoCut = 0.1;                // 10%
};

// destructor
AnalysisDD4hep::~AnalysisDD4hep() {
  delete AN;
};

void AnalysisDD4hep::process_event()
{

  PIDtoEnum_ = AN->GetPIDMap();
      
  // instantiate objects
  kin = new Kinematics(eleBeamEn,ionBeamEn,crossingAngle);
  kinTrue = new Kinematics(eleBeamEn, ionBeamEn, crossingAngle);
  ST = new SimpleTree("tree",kin);
  weight = new WeightsUniform();
  weightJet = new WeightsUniform();

  TChain *chain = new TChain("events");
  for(int i=0; i<(int)infiles.size(); i++)
    {
      cout << infiles[i].Data() << endl;
      chain->Add(infiles[i].Data());
    }

  // FIXME: replace it with ExRootTreeReader::UseBranch()?
  TTreeReader tr(chain);

  // Truth
  TTreeReaderArray<Int_t>    mcparticles2_pdgID(tr,     "mcparticles2.pdgID");
  TTreeReaderArray<Double_t> mcparticles2_psx(tr,       "mcparticles2.ps.x");
  TTreeReaderArray<Double_t> mcparticles2_psy(tr,       "mcparticles2.ps.y");
  TTreeReaderArray<Double_t> mcparticles2_psz(tr,       "mcparticles2.ps.z");
  TTreeReaderArray<Int_t>    mcparticles2_status(tr,    "mcparticles2.status");
  TTreeReaderArray<Int_t>    mcparticles2_genStatus(tr, "mcparticles2.genStatus");
  TTreeReaderArray<Double_t> mcparticles2_mass(tr,      "mcparticles2.mass");

  // Reco
  TTreeReaderArray<int> ReconstructedParticles_pid(tr, "ReconstructedParticles.pid");
  TTreeReaderArray<float> ReconstructedParticles_energy(tr, "ReconstructedParticles.energy");
  TTreeReaderArray<float> ReconstructedParticles_p_x(tr, "ReconstructedParticles.p.x");
  TTreeReaderArray<float> ReconstructedParticles_p_y(tr, "ReconstructedParticles.p.y");
  TTreeReaderArray<float> ReconstructedParticles_p_z(tr, "ReconstructedParticles.p.z");
  TTreeReaderArray<float> ReconstructedParticles_mass(tr, "ReconstructedParticles.mass");

  //HcalEndcap
  TTreeReaderArray<float>  HcalEndcapPClusters_energy(tr, "HcalEndcapPClusters.energy");
  TTreeReaderArray<float> HcalEndcapPClusters_x(tr,      "HcalEndcapPClusters.position.x");
  TTreeReaderArray<float> HcalEndcapPClusters_y(tr,      "HcalEndcapPClusters.position.y");
  TTreeReaderArray<float> HcalEndcapPClusters_z(tr,      "HcalEndcapPClusters.position.z");
  TTreeReaderArray<float>  HcalEndcapPClusters_theta(tr,  "HcalEndcapPClustersInfo.polar.theta");
  TTreeReaderArray<float>  HcalEndcapPClusters_phi(tr,    "HcalEndcapPClustersInfo.polar.phi");

  //HcalEndcap
  TTreeReaderArray<float> HcalEndcapNClusters_energy(tr, "HcalEndcapNClusters.energy");
  TTreeReaderArray<float> HcalEndcapNClusters_x(tr,      "HcalEndcapNClusters.position.x");
  TTreeReaderArray<float> HcalEndcapNClusters_y(tr,      "HcalEndcapNClusters.position.y");
  TTreeReaderArray<float> HcalEndcapNClusters_z(tr,      "HcalEndcapNClusters.position.z");
  TTreeReaderArray<float> HcalEndcapNClusters_theta(tr,  "HcalEndcapNClustersInfo.polar.theta");
  TTreeReaderArray<float> HcalEndcapNClusters_phi(tr,    "HcalEndcapNClustersInfo.polar.phi");

  //HcalBarrel
  TTreeReaderArray<float> HcalBarrelClusters_energy(tr, "HcalBarrelClusters.energy");
  TTreeReaderArray<float> HcalBarrelClusters_x(tr,      "HcalBarrelClusters.position.x");
  TTreeReaderArray<float> HcalBarrelClusters_y(tr,      "HcalBarrelClusters.position.y");
  TTreeReaderArray<float> HcalBarrelClusters_z(tr,      "HcalBarrelClusters.position.z");
  TTreeReaderArray<float> HcalBarrelClusters_theta(tr,  "HcalBarrelClustersInfo.polar.theta");
  TTreeReaderArray<float> HcalBarrelClusters_phi(tr,    "HcalBarrelClustersInfo.polar.phi");

  //Ecal
  TTreeReaderArray<float> EcalEndcapPClusters_energy(tr, "EcalEndcapPClusters.energy");
  TTreeReaderArray<float> EcalEndcapPClusters_x(tr,      "EcalEndcapPClusters.position.x");
  TTreeReaderArray<float> EcalEndcapPClusters_y(tr,      "EcalEndcapPClusters.position.y");
  TTreeReaderArray<float> EcalEndcapPClusters_z(tr,      "EcalEndcapPClusters.position.z");
  TTreeReaderArray<float> EcalEndcapPClusters_theta(tr,  "EcalEndcapPClustersInfo.polar.theta");
  TTreeReaderArray<float> EcalEndcapPClusters_phi(tr,    "EcalEndcapPClustersInfo.polar.phi");

  TTreeReaderArray<float> EcalEndcapNClusters_energy(tr, "EcalEndcapNClusters.energy");
  TTreeReaderArray<float> EcalEndcapNClusters_x(tr,      "EcalEndcapNClusters.position.x");
  TTreeReaderArray<float> EcalEndcapNClusters_y(tr,      "EcalEndcapNClusters.position.y");
  TTreeReaderArray<float> EcalEndcapNClusters_z(tr,      "EcalEndcapNClusters.position.z");
  //watch out the branch name
  TTreeReaderArray<float>  EcalEndcapNClusters_theta(tr,  "EcalEndcapNClusterInfo.polar.theta");
  TTreeReaderArray<float>  EcalEndcapNClusters_phi(tr,    "EcalEndcapNClusterInfo.polar.phi");

  TTreeReaderArray<float> EcalBarrelClusters_energy(tr, "EcalBarrelImagingClusters.energy");
  TTreeReaderArray<float> EcalBarrelClusters_x(tr,      "EcalBarrelImagingClusters.position.x");
  TTreeReaderArray<float> EcalBarrelClusters_y(tr,      "EcalBarrelImagingClusters.position.y");
  TTreeReaderArray<float> EcalBarrelClusters_z(tr,      "EcalBarrelImagingClusters.position.z");
  TTreeReaderArray<float> EcalBarrelClusters_theta(tr,  "EcalBarrelImagingClustersInfo.polar.theta");
  TTreeReaderArray<float>  EcalBarrelClusters_phi(tr,    "EcalBarrelImagingClustersInfo.polar.phi");

  TTreeReader::EEntryStatus entrystats = tr.SetEntry(0);

  // open output file
  // FIXME: replace name variable
  outFile = new TFile(outfilePrefix,"RECREATE");

  // number of bins

  const Int_t NptBins = AN->BinScheme("pt")->GetNumBins();
  const Int_t NxBins = AN->BinScheme("x")->GetNumBins();
  const Int_t NzBins = AN->BinScheme("z")->GetNumBins();
  const Int_t NqBins = AN->BinScheme("q2")->GetNumBins();
  const Int_t NyBins = AN->BinScheme("y")->GetNumBins();
  const Int_t NfinalStateBins = AN->BinScheme("finalState")->GetNumBins();
  const Int_t NptjetBins = AN->BinScheme("pt_jet")->GetNumBins();
  const Int_t NzjetBins = AN->BinScheme("z_jet")->GetNumBins();
  const Int_t NrecMethodBins = AN->BinScheme("recMethod")->GetNumBins();

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
          if(AN->CheckDiagonal(bpt,bx,bz,bq)) continue; // diagonalizer
          for(int by=0; by<NyBins; by++) { // - loop over y bins
            for(int bfs=0; bfs<NfinalStateBins; bfs++) { // - loop over final states

              // set Histos name and title
              histosN = AN->GetHistosName (bpt,bx,bz,bq,by,bfs);
              histosT = AN->GetHistosTitle(bpt,bx,bz,bq,by,bfs);

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
              HS->DefineHist1D("xRes","x - x_{true}","", NBINS, -2, 2);
              HS->DefineHist1D("yRes","y - y_{true}","", NBINS, -2, 2);
              HS->DefineHist1D("Q2Res","Q2 - Q2_{true}","", NBINS, -2, 2);
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
              HS->AddCutDef(AN->BinScheme("pt")->Cut(bpt));
              HS->AddCutDef(AN->BinScheme("x")->Cut(bx));
              HS->AddCutDef(AN->BinScheme("z")->Cut(bz));
              HS->AddCutDef(AN->BinScheme("q2")->Cut(bq));
              HS->AddCutDef(AN->BinScheme("y")->Cut(by));
              HS->AddCutDef(AN->BinScheme("finalState")->Cut(bfs));

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

            histosN = AN->GetHistosNameJets(bpt, bz, bx, bq, by);
            histosT = AN->GetHistosTitleJets(bpt, bz, bx, bq, by);

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
            HS->AddCutDef(AN->BinScheme("pt_jet")->Cut(bpt));
            HS->AddCutDef(AN->BinScheme("z_jet")->Cut(bz));
            HS->AddCutDef(AN->BinScheme("x")->Cut(bx));
            HS->AddCutDef(AN->BinScheme("q2")->Cut(bq));
            HS->AddCutDef(AN->BinScheme("y")->Cut(by));
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
              histosN = AN->GetHistosNameBreitJets(bpt, bz, bx, bq, by, brec);
              histosT = AN->GetHistosTitleBreitJets(bpt, bz, bx, bq, by, brec);
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
              HS->AddCutDef(AN->BinScheme("pt_jet")->Cut(bpt));
              HS->AddCutDef(AN->BinScheme("z_jet")->Cut(bz));
              HS->AddCutDef(AN->BinScheme("x")->Cut(bx));
              HS->AddCutDef(AN->BinScheme("q2")->Cut(bq));
              HS->AddCutDef(AN->BinScheme("y")->Cut(by));
              HS->AddCutDef(AN->BinScheme("recMethod")->Cut(brec));
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

  int noele = 0;
  // event loop =========================================================
  cout << "begin event loop..." << endl;
  int nevt = 0;
  while(tr.Next())
    {
      if(nevt%100000==0) cout << nevt << " events..." << endl;
      nevt++;      

      std::vector<Particles> mcpart;
      double maxP = 0;
      for(int imc=0; imc<mcparticles2_pdgID.GetSize(); imc++)
	{
	  int pid_ = mcparticles2_pdgID[imc];
	  double px_ = mcparticles2_psx[imc];
	  double py_ = mcparticles2_psy[imc];
	  double pz_ = mcparticles2_psz[imc];
	  double mass_ = mcparticles2_mass[imc]; // in GeV
	  double p_ = sqrt(pow(mcparticles2_psx[imc],2) + pow(mcparticles2_psy[imc],2) + pow(mcparticles2_psz[imc],2));

	  // genStatus 4: beam particle 1: final state 
	  if(mcparticles2_genStatus[imc] == 1)
	    {
	      Particles part;
	      part.pid = pid_;
	      part.vecPart.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
	      mcpart.push_back(part);

	      if(mcparticles2_pdgID[imc] == 11)
		{
		  if(p_ > maxP)
		    {
		      maxP = p_;
		      kinTrue->vecElectron.SetPxPyPzE(mcparticles2_psx[imc],
						      mcparticles2_psy[imc],
						      mcparticles2_psz[imc],
						      sqrt(p_*p_ + mass_*mass_));
		    }
		}// if electron
	    }//
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
      // HcalEndcapP
      for(int icl=0; icl<HcalEndcapPClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = HcalEndcapPClusters_energy[icl]; // use edep?
	  clus->x = HcalEndcapPClusters_x[icl]; // use edep?
	  clus->y = HcalEndcapPClusters_y[icl]; // use edep?
	  clus->z = HcalEndcapPClusters_z[icl]; // use edep?
	  clus->theta = HcalEndcapPClusters_theta[icl]; // use edep?
	  clus->phi = HcalEndcapPClusters_phi[icl]; // use edep?
	  v_hcal_clusters.push_back(clus);
	}

      // HcalEndcapN
      for(int icl=0; icl<HcalEndcapNClusters_x.GetSize(); icl++)
	{
	  Clusters* clus = new Clusters();
	  clus->E = HcalEndcapNClusters_energy[icl]; // use edep?
	  clus->x = HcalEndcapNClusters_x[icl]; // use edep?
	  clus->y = HcalEndcapNClusters_y[icl]; // use edep?
	  clus->z = HcalEndcapNClusters_z[icl]; // use edep?
	  clus->theta = HcalEndcapNClusters_theta[icl]; // use edep?
	  clus->phi = HcalEndcapNClusters_phi[icl]; // use edep?
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

      // calculate true DIS kinematics
      kinTrue->CalculateDISbyElectron(); // generated (truth)

      // find scattered electron
      int electron_index = find_electron(v_ecal_clusters, v_hcal_clusters, fEThreshold);
      if(electron_index < 0)
	{
	  //cout << nevt << " No scattered electron found.. skip this event" << electron_index << endl;
	  //	  cout << kinTrue->x << " " << kinTrue->Q2 << " " << kinTrue->y << " " << kinTrue->vecElectron.Pt() << " " << kinTrue->vecElectron.P() << " " << kinTrue->vecElectron.Eta() << endl;
	  noele++;
	  continue;
	}

      double electron_E     = v_ecal_clusters[electron_index]->E;
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

      TLorentzVector v_had;
      double hpx=0; 
      double hpy=0;
      for(int itrk=0; itrk<(int)ReconstructedParticles_pid.GetSize(); itrk++)
	{
	  // FIXME: pid is using the true information
	  // Add PID smearing
	  int pid_ = ReconstructedParticles_pid[itrk];

	  // pid==0: reconstructed tracks with no matching truth pid
	  if(pid_ == 0) continue;

	  auto kv = PIDtoEnum_.find(pid_);
	  if(kv!=PIDtoEnum_.end()) bFinalState = kv->second;
	  else continue;

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

	  // find the true info
	  double mineta = 4.0;
	  for(int imc=0; imc<(int)mcpart.size(); imc++)
	    {
	      if(pid_ == mcpart[imc].pid)
		{
		  double deta = abs(kin->vecHadron.Eta() - mcpart[imc].vecPart.Eta());
		  if( deta < mineta )
		    {
		      mineta = deta;
		      kinTrue->vecHadron = mcpart[imc].vecPart;
		    }
		}
	    }

	  kinTrue->CalculateHadronKinematics();

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
	    CheckBins( AN->BinScheme("pt"), v_pt, kin->pT );
	    CheckBins( AN->BinScheme("x"),  v_x,  kin->x );
	    CheckBins( AN->BinScheme("z"),  v_z,  kin->z );
	    CheckBins( AN->BinScheme("q2"), v_q,  kin->Q2 );
	    CheckBins( AN->BinScheme("y"),  v_y,  kin->y );

	    // build list of histogram sets to fill
	    histSetFillList.clear();
	    for(int bpt : v_pt) {
	      for(int bx : v_x) {
		for(int bz : v_z) {
		  for(int bq : v_q) {
		    for(int by : v_y) {
		      if(!AN->CheckDiagonal(bpt,bx,bz,bq)) {
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
	      H->Hist("xRes")->Fill(kin->x - kinTrue->x,w);
	      H->Hist("yRes")->Fill(kin->y - kinTrue->y,w);
	      H->Hist("Q2Res")->Fill(kin->Q2 - kinTrue->Q2,w);
	    };


	    // fill simple tree (not binned)
	    // TODO: consider adding a `finalState` cut
	    if( AN->writeSimpleTree && histSetFillList.size()>0 ) ST->FillTree(w);
	
	  }//if cut
	}//trk loop

      //Hadronic reconstruction method
      kin->sigmah = (kin->vecElectron.E() - kin->vecElectron.Pz());
      kin->Pxh = hpx - kin->vecElectron.Px();
      kin->Pyh = hpy - kin->vecElectron.Py();

    }// tree reader loop

  cout << "Total no scattered electron found: " << noele << endl;
  cout << "end event loop" << endl;
  // event loop end =========================================================                 
  // calculate integrated luminosity  
  Double_t lumi = wTotal/xsecTot; // [nb^-1]                                                      

  cout << "Integrated Luminosity:       " << lumi << "/nb" << endl;
  cout << sep << endl;

  // calculate cross sections 
  // TODO: generalize (`if (name contains "xsec") ...`) 
  for(Histos *H : histSetList) {
    H->Hist("Q_xsec")->Scale(1./lumi);
  };
  // print yields in each bin
  cout << sep << endl << "Histogram Entries:" << endl;
  for(Histos *H : histSetList) {
    cout << H->GetSetTitle() << " ::: "
	 << H->Hist("Q2vsX")->GetEntries()
	 << endl;
  };
  
  // write histograms 
  cout << sep << endl;
  outFile->cd();
  if(AN->writeSimpleTree) ST->WriteTree();
  for(Histos *H : histSetList) H->WriteHists(outFile);
  for(Histos *H : histSetList) H->Write();
  for(Histos *H : histSetListJets) H->WriteHists(outFile);
  for(Histos *H : histSetListJets) H->Write();
#if INCCENTAURO == 1
  for(Histos *H : histSetListBreitJets) H->WriteHists(outfile);
  for(Histos *H : histSetListBreitJets) H->Write();
#endif
  // write binning schemes
  binSchemes_ = AN->GetBinSchemes();
  for(auto const &kv : binSchemes_) kv.second->Write(kv.first+"_bins");
  // close output
  outFile->Close();
  cout << outfileName << " written." << endl;
  

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
      if( dphi > TMath::Pi() ) dphi = dphi - 2*TMath::Pi();
      double dr = sqrt( pow(dphi,2) + pow((clus_eta - cone_eta), 2) );

      // get E_cone
      if(dr < fIsoR)
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
      double clus_E = icluster.E;
      double clus_theta = icluster.theta;
      double clus_phi   = icluster.phi;
      double clus_eta   = -log(tan(clus_theta/2.0));
      double clus_pt    = clus_E * sin(clus_theta);

      if(clus_E < e_threshold)
	continue;
      if(icluster.theta * 180./TMath::Pi() < 2)
	continue;

      // FIXME hard-coded threshold (currently not used)
      double clus_ecal_iso = isolation(clus_theta, clus_phi, ecal_cluster, 0.1);
      double clus_hcal_iso = isolation(clus_theta, clus_phi, hcal_cluster, 0.1);
      double clus_iso = clus_ecal_iso + clus_hcal_iso - clus_E; // subtract the electron energy

      if(clus_iso > clus_E*fIsoCut)
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

void AnalysisDD4hep::CheckBins(BinSet *bs, std::vector<int> &v, Double_t var) {
  v.clear();
  for(int b=0; b<bs->GetNumBins(); b++) {
    if(bs->Cut(b)->CheckCut(var)) v.push_back(b);
  };
};
