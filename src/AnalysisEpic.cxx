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
  for(unsigned e=0; e<podioReader.getEntries(); e++) {
    if(e%10000==0) cout << e << " events..." << endl;
    cout << endl << "EVENT " << e << " ===============================================" << endl;

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();
    double mcpartElectronP   = 0.0;
    bool double_counted_beam = false;
    int num_ele_beams    = 0;
    int num_ion_beams    = 0;
    int num_mc_electrons = 0;
    
    // read particle collections for this event
    auto& mcparts  = evStore.get<edm4hep::MCParticleCollection>("MCParticles");
    auto& recparts = evStore.get<edm4hep::ReconstructedParticleCollection>("ReconstructedParticles");

    // data objects
    edm4hep::MCParticle mcpartEleBeam;
    edm4hep::MCParticle mcpartIonBeam;
    edm4hep::MCParticle mcpartElectron;

    // loop over generated particles
    cout << endl << "MCParticles: ---------------------------------------" << endl;
    for(auto mcpart : mcparts) {

      // print out this MCParticle
      // PrintParticle(mcpart);

      // filter for beam particles
      if(mcpart.getGeneratorStatus() == constants::statusBeam) {
        switch(mcpart.getPDG()) {
          case constants::pdgElectron:
            if(num_ele_beams>0) double_counted_beam = true;
            mcpartEleBeam = mcpart;
            num_ele_beams++;
            break;
          case constants::pdgProton:
            if(num_ion_beams>0) double_counted_beam = true;
            mcpartIonBeam = mcpart;
            num_ion_beams++;
            break;
          default:
            ErrorPrint(Form("WARNING: Unknown beam particle with PDG=%d",mcpart.getPDG()));
        }
      }

      // filter for scattered electron: select the one with the highest |p|
      if(mcpart.getGeneratorStatus() == constants::statusFinal) {
        if(mcpart.getPDG() == constants::pdgElectron) {
          auto eleP = edm4hep::utils::p(mcpart);
          if(eleP>mcpartElectronP) {
            mcpartElectron  = mcpart;
            mcpartElectronP = eleP;
            num_mc_electrons++;
          }
        }
      }

    } // loop over generated particles

    // check for found generated particles
    if(num_ele_beams==0)    { ErrorPrint("WARNING: missing MC electron beam");      continue; };
    if(num_ion_beams==0)    { ErrorPrint("WARNING: missing MC ion beam");           continue; };
    if(num_mc_electrons==0) { ErrorPrint("WARNING: missing scattered electron");    continue; };
    if(double_counted_beam) { ErrorPrint("WARNING: found multiple beam particles"); continue; };

    // print beam particles
    cout << endl << "GENERATED BEAMS ------------------------------------" << endl;
    PrintParticle(mcpartEleBeam);
    PrintParticle(mcpartIonBeam);
    cout << endl << "GENERATED SCATTERED ELECTRON -----------------------" << endl;
    PrintParticle(mcpartElectron);

    // loop over reconstructed particles
    cout << endl << "ReconstructedParticles: ----------------------------" << endl;
    for(const auto& recpart : recparts) {

      // print out this ReconstructedParticle
      PrintParticle(recpart);

      // for(const auto& track : P.getTracks()) {
      //   // ...
      // }
      // for(const auto& cluster : P.getClusters()) {
      //   // ...
      // }

    } // loop over reconstructed particles


    // read kinematics calculations from upstream /////////////////////////
    // TODO: cross check these with our calculations from `Kinematics`
    cout << endl << "KINEMATICS, calculated from upstream:" << endl;
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
  cout << "  Momentum: ("
    << P.getMomentum().x << ", "
    << P.getMomentum().y << ", "
    << P.getMomentum().z << ")" << endl;
  cout << "  p=|Momentum|: " << edm4hep::utils::p(P) << endl;
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
  cout << "  Momentum: ("
    << P.getMomentum().x << ", "
    << P.getMomentum().y << ", "
    << P.getMomentum().z << ")" << endl;
  cout << "  p=|Momentum|: " << edm4hep::utils::p(P) << endl;
  cout << "  Energy:       " << P.getEnergy() << endl;
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
