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
};


void AnalysisDD4hep::process_event()
{
  // setup
  Prepare();

  // read dd4hep tree
  TChain *chain = new TChain("events");
  for(TString in : infiles) chain->Add(in);

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


  int noele = 0;
  // event loop =========================================================
  cout << "begin event loop..." << endl;
  Long64_t nevt = 0;
  while(tr.Next())
    {
      if(nevt%10000==0) cout << nevt << " events..." << endl;
      nevt++;      
      if(nevt>maxEvents) break;

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

      // calculate true DIS kinematics
      kinTrue->CalculateDIS(reconMethod); // generated (truth)

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

      // find scattered electron
      int electron_index = find_electron(v_ecal_clusters, v_hcal_clusters, fEThreshold);
      if(electron_index < 0)
	{
	  noele++;
	  continue;
	}

      // Set electron kinematics
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

      // track loop
      std::vector<Particles> recopart;
      double hpx=0; 
      double hpy=0;
      for(int itrk=0; itrk<(int)ReconstructedParticles_pid.GetSize(); itrk++)
	{
	  // FIXME: pid is using the true information
	  // Add PID smearing
	  int pid_ = ReconstructedParticles_pid[itrk];

	  // pid==0: reconstructed tracks with no matching truth pid
	  if(pid_ == 0) continue;

	  Particles part;
	  part.pid = pid_;

	  TLorentzVector v_temp;
	  double reco_E = ReconstructedParticles_energy[itrk];
	  double reco_px = ReconstructedParticles_p_x[itrk];
	  double reco_py = ReconstructedParticles_p_y[itrk];
	  double reco_pz = ReconstructedParticles_p_z[itrk];
	  double reco_mass = ReconstructedParticles_mass[itrk];
	  double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
	  v_temp.SetPxPyPzE(reco_px, reco_py, reco_pz, reco_E);

	  part.vecPart.SetPxPyPzE(reco_px, 
				  reco_py, 
				  reco_pz, 
				  sqrt(reco_p*reco_p + reco_mass*reco_mass));

	  recopart.push_back(part);
	  hpx += reco_px;
	  hpy += reco_py;
	}

      //Hadronic reconstruction 
      kin->sigmah = (kin->vecElectron.E() - kin->vecElectron.Pz());
      kin->Pxh = hpx - kin->vecElectron.Px();
      kin->Pyh = hpy - kin->vecElectron.Py();

      // calculate DIS kinematics
      kin->CalculateDIS(reconMethod); // reconstructed

      for(auto trk : recopart)
	{
	  int pid_ = trk.pid;

    // final state cut
    // - check PID, to see if it's a final state we're interested in for
    //   histograms; if not, proceed to next track
    auto kv = PIDtoFinalState.find(pid_);
    if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
    if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

	  kin->vecHadron = trk.vecPart;

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

    // asymmetry injection
    //kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
    //kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
    //kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

	  wTrack = weight->GetWeight(*kinTrue);
	  wTrackTotal += wTrack;

	  // apply cuts
	  if(kin->CutFull()) {

      // fill track histograms in activated bins
      FillHistosTracks();

      // fill simple tree
      // - not binned
      // - `activeEvent` is only true if at least one bin gets filled for this track
      // - TODO [critical]: add a `finalState` cut (also needed in AnalysisDelphes)
      if( writeSimpleTree && activeEvent ) ST->FillTree(wTrack);

	  }//if cut
	}//trk loop

    }// tree reader loop

  cout << "Total no scattered electron found: " << noele << endl;
  cout << "end event loop" << endl;
  // event loop end =========================================================                 
  
  // finish execution
  Finish();

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

