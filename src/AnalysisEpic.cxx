#include "AnalysisEpic.h"

using std::cout;
using std::cerr;
using std::endl;

AnalysisEpic::AnalysisEpic(
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
      ) 
{};

AnalysisEpic::~AnalysisEpic() {};

void AnalysisEpic::Execute()
{
  // setup
  Prepare();

  // produce flat list of files from `infiles`
  std::vector<std::string> infilesFlat;
  for(const auto fileList : infiles)
    for(const auto fileName : fileList)
      infilesFlat.push_back(fileName);

  // create PODIO event store
  podioReader.openFiles(infilesFlat);
  evStore.setReader(&podioReader);
  ENT = podioReader.getEntries();
  if(maxEvents>0) ENT = std::min(maxEvents,ENT);
  
  // calculate Q2 weights
  CalculateEventQ2Weights();

  // list of reconstruction methods produced upstream
  const std::vector<std::string> upstreamReconMethodList = {
    "Truth",
    "Electron",
    "DA",
    "JB",
    "Sigma"
  };

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  for(unsigned e=0; e<ENT; e++) {
    if(e%10000==0) cout << e << " events..." << endl;
    if(verbose)
      PrintHeader(Form("EVENT %d ===============================================",e));

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();
    double mcPartElectronP   = 0.0;
    bool double_counted_beam = false;
    int num_ele_beams    = 0;
    int num_ion_beams    = 0;
    int num_sim_electrons = 0;
    int num_rec_electrons = 0;
    
    // read particle collections for this event
    auto& simParts  = evStore.get<edm4hep::MCParticleCollection>("MCParticles");
    // auto& recParts = evStore.get<edm4hep::ReconstructedParticleCollection>("ReconstructedParticles");
    auto& mcRecLinks = evStore.get<edm4hep::MCRecoParticleAssociation>("ReconstructedParticlesAssoc");

    // data objects
    edm4hep::MCParticle mcPartEleBeam;
    edm4hep::MCParticle mcPartIonBeam;
    edm4hep::MCParticle mcPartElectron;
    std::set<edm4hep::ReconstructedParticle> recPartsToAnalyze;

    // loop over generated particles
    PrintHeader("MCParticles: ---------------------------------------");
    for(auto simPart : simParts) {

      // print out this MCParticle
      // if(verbose) PrintParticle(simPart);

      // generated particle properties
      auto pid = simPart.getPDG();

      // add to Hadronic Final State (HFS) sums
      kinTrue->AddToHFS(GetP4(simPart));

      // filter for beam particles
      if(simPart.getGeneratorStatus() == constants::statusBeam) {
        switch(pid) {
          case constants::pdgElectron:
            if(num_ele_beams>0) double_counted_beam = true;
            mcPartEleBeam = simPart;
            num_ele_beams++;
            break;
          case constants::pdgProton:
            if(num_ion_beams>0) double_counted_beam = true;
            mcPartIonBeam = simPart;
            num_ion_beams++;
            break;
          default:
            ErrorPrint(Form("WARNING: Unknown beam particle with PDG=%d",pid));
        }
      }

      // filter for scattered electron: select the one with the highest |p|
      if(simPart.getGeneratorStatus() == constants::statusFinal) {
        if(pid == constants::pdgElectron) {
          auto eleP = edm4hep::utils::p(simPart);
          if(eleP>mcPartElectronP) {
            mcPartElectron  = simPart;
            mcPartElectronP = eleP;
            num_sim_electrons++;
          }
        }
      }

    } // loop over generated particles

    // check for found generated particles
    if(num_ele_beams==0)     { ErrorPrint("WARNING: missing MC electron beam");      continue; };
    if(num_ion_beams==0)     { ErrorPrint("WARNING: missing MC ion beam");           continue; };
    if(num_sim_electrons==0) { ErrorPrint("WARNING: missing scattered electron");    continue; };
    if(double_counted_beam)  { ErrorPrint("WARNING: found multiple beam particles"); continue; };

    // set Kinematics 4-momenta
    kinTrue->vecEleBeam  = GetP4(mcPartEleBeam);
    kinTrue->vecIonBeam  = GetP4(mcPartIonBeam);
    kinTrue->vecElectron = GetP4(mcPartElectron);

    // print beam particles
    if(verbose) {
      PrintHeader("GENERATED BEAMS ------------------------------------");
      PrintParticle(mcPartEleBeam);
      PrintParticle(mcPartIonBeam);
      PrintHeader("GENERATED SCATTERED ELECTRON -----------------------");
      PrintParticle(mcPartElectron);
    }

    // loop over reconstructed particles
    /*
    if(verbose) PrintHeader("ReconstructedParticles: ----------------------------");
    for(const auto& recPart : recParts) {

      // print out this ReconstructedParticle
      if(verbose) PrintParticle(recPart);

      // for(const auto& track : recPart.getTracks()) {
      //   // ...
      // }
      // for(const auto& cluster : recPart.getClusters()) {
      //   // ...
      // }

    } // loop over reconstructed particles
    */


    // loop over associations: MC particle <-> Reconstructed particle
    for(const auto& link : mcRecLinks ) {
      auto simPart = link.getSim();
      auto recPart = link.getRec();

      // get PID
      int pid = 0;
      // get reconstructed PID
      if(recPart.getParticleIDUsed().isAvailable()) // FIXME: is always false, not yet available upstream
        pid = P.getParticleIDUsed().getPDG();
      if(pid==0) {
        // continue; // realistic: skip this particle
        pid = simPart.getPDG(); // unrealistic fix: if reconstructed PID is unavailable, use MC PID
      }

      // add to list of reconstructed particles to analyze, and to the HFS
      recPartsToAnalyze.insert(recPart);
      kin->AddToHFS(GetP4(recPart));

      // find scattered electron, by matching to truth
      // FIXME: not realistic
      // FIXME: does `simP==mcPartElectron` actually work !?
      if(pid==constants::pdgElectron && simP==mcPartElectron) {
        num_rec_electrons++;
        kin->vecElectron = GetP4(recPart);
      }
    }

    // check for found reconstructed particles
    if(num_rec_electrons == 0) { ErrorPrint("WARNING: reconstructed scattered electron not found"); continue; };
    if(num_rec_electrons >  1) { ErrorPrint("WARNING: found more than 1 reconstructed scattered electron"); };

    // subtract electron from hadronic final state variables
    kin->SubtractElectronFromHFS();
    kinTrue->SubtractElectronFromHFS();

    // skip the event if there are no reconstructed particles (other than the
    // electron), otherwise hadronic recon methods will fail
    if(kin->countHadrons == 0) { ErrorPrint("WARNING: no hadrons"); };
    
    // calculate DIS kinematics
    if(!(kin->CalculateDIS(reconMethod))) continue; // reconstructed
    if(!(kinTrue->CalculateDIS(reconMethod))) continue; // generated (truth)


    //
    //
    // TODO: stopped here
    //
    //


    // read kinematics calculations from upstream /////////////////////////
    // TODO: cross check these with our calculations from `Kinematics`
    PrintHeader("KINEMATICS, calculated from upstream:");
    printf("  %10s %8s %8s %8s %8s %8s %8s\n", "", "x", "Q2", "W", "y", "nu", "elec?");
    for(const auto upstreamReconMethod : upstreamReconMethodList)
      for(const auto& calc : evStore.get<eicd::InclusiveKinematicsCollection>("InclusiveKinematics"+upstreamReconMethod) )
        printf("  %10s %8.5f %8.2f %8.2f %8.5f %8.2f %8d\n",
            upstreamReconMethod.c_str(),
            calc.getX(),
            calc.getQ2(),
            calc.getW(),
            calc.getY(),
            calc.getNu(),
            calc.getScat().isAvailable()? 1:0
            );


    // next event //////////////////////////////////
    evStore.clear();
    podioReader.endOfEvent();

  } // event loop
  cout << "end event loop" << endl;

  // finish execution
  podioReader.closeFile();
  Finish();
}


// particle printers //////////////////////////////////////////////

void AnalysisEpic::PrintParticle(const edm4hep::MCParticle& P) { 
  cout << endl;
  cout << "  PDG: " << P.getPDG() << endl;
  cout << "  Status: " << P.getGeneratorStatus() << endl;
  cout << "  Vertex:   ("
    << P.getVertex().x << ", "
    << P.getVertex().y << ", "
    << P.getVertex().z << ")" << endl;
  cout << "  p=|Momentum|: " << edm4hep::utils::p(P) << endl;
  cout << "  Energy:       " << P.getEnergy() << endl;
  cout << "  3-Momentum: ("
    << P.getMomentum().x << ", "
    << P.getMomentum().y << ", "
    << P.getMomentum().z << ")" << endl;
  cout << "  4-Momentum: ";
  GetP4(P).Print();
  cout << "  pT_lab:       " << edm4hep::utils::pT(P) << endl;
  cout << "  Parents:" << endl;
  for(const auto& parent : P.getParents())
    cout << "    PDG: " << parent.getPDG() << endl;
  cout << "  Daughters:" << endl;
  for(const auto& daughter : P.getDaughters())
    cout << "    PDG: " << daughter.getPDG() << endl;
}

void AnalysisEpic::PrintParticle(const edm4hep::ReconstructedParticle& P) {
  cout << endl;
  if(P.getParticleIDUsed().isAvailable())
    cout << "  PDG: " << P.getParticleIDUsed().getPDG() << endl;
  else cout << "  PDG: ???" << endl;
  // if(P.getStartVertex().isAvailable()) // FIXME: sometimes segfaults?
  //   cout << "  StartVertex: ("
  //     << P.getStartVertex().getPosition().x << ", "
  //     << P.getStartVertex().getPosition().y << ", "
  //     << P.getStartVertex().getPosition().z << ")" << endl;
  // else cout << "  StartVertex: ???" << endl;
  cout << "  p=|Momentum|: " << edm4hep::utils::p(P) << endl;
  cout << "  Energy:       " << P.getEnergy() << endl;
  cout << "  3-Momentum: ("
    << P.getMomentum().x << ", "
    << P.getMomentum().y << ", "
    << P.getMomentum().z << ")" << endl;
  cout << "  4-Momentum: ";
  GetP4(P).Print();
  cout << "  Mass:         " << P.getMass() << endl;
  cout << "  Charge:       " << P.getCharge() << endl;
  cout << "  # of clusters: "   << P.clusters_size()    << endl;
  cout << "  # of tracks:   "   << P.tracks_size()      << endl;
  cout << "  # of PIDs:     "   << P.particleIDs_size() << endl;
  cout << "  # of combined reconstructed parts: " << P.particles_size() << endl;
  // for(const auto& track : P.getTracks()) {
  //   // ...
  // }
  // for(const auto& cluster : P.getClusters()) {
  //   // ...
  // }
}

void AnalysisEpic::PrintHeader(TString msg) {
  cout << endl << msg << endl;
};
