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
  for(auto fileList : infiles)
    for(auto fileName : fileList)
      infilesFlat.push_back(fileName);

  // create PODIO event store
  podioReader.openFiles(infilesFlat);
  evStore.setReader(&podioReader);
  
  // calculate Q2 weights
  CalculateEventQ2Weights();

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  for(unsigned e=0; e<podioReader.getEntries(); e++) {
    if(e%10000==0) cout << e << " events..." << endl;
    cout << endl << "EVENT " << e << " ===============================================" << endl;

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();
    
    // particle collections
    auto& mcparts  = evStore.get<edm4hep::MCParticleCollection>("MCParticles");
    auto& recparts = evStore.get<edm4hep::ReconstructedParticleCollection>("ReconstructedParticles");

    // generated particles loop
    cout << endl << "MCParticles: ------------" << endl;
    for(const auto& mcpart : mcparts) {
      cout << endl;
      cout << "  PDG: " << mcpart.getPDG() << endl;
      cout << "  Status: " << mcpart.getGeneratorStatus() << endl;
      cout << "  Vertex:   ("
        << mcpart.getVertex().x << ", "
        << mcpart.getVertex().y << ", "
        << mcpart.getVertex().z << ")" << endl;
      cout << "  Momentum: ("
        << mcpart.getMomentum().x << ", "
        << mcpart.getMomentum().y << ", "
        << mcpart.getMomentum().z << ")" << endl;
      cout << "  p=|Momentum|: " << edm4hep::utils::p(mcpart) << endl;
      cout << "  pT_lab:       " << edm4hep::utils::pT(mcpart) << endl;
      cout << "  Parents:" << endl;
      for(const auto& parent : mcpart.getParents())
        cout << "    PDG: " << parent.getPDG() << endl;
      cout << "  Daughters:" << endl;
      for(const auto& daughter : mcpart.getDaughters())
        cout << "    PDG: " << daughter.getPDG() << endl;
    }

    // reconstructed particles loop
    cout << endl << "ReconstructedParticles: ------------" << endl;
    for(const auto& recpart : recparts) {
      cout << endl;
      if(recpart.getParticleIDUsed().isAvailable())
        cout << "  PDG: " << recpart.getParticleIDUsed().getPDG() << endl;
      else cout << "  PDG: ???" << endl;
      // if(recpart.getStartVertex().isAvailable()) // FIXME: sometimes segfaults?
      //   cout << "  StartVertex: ("
      //     << recpart.getStartVertex().getPosition().x << ", "
      //     << recpart.getStartVertex().getPosition().y << ", "
      //     << recpart.getStartVertex().getPosition().z << ")" << endl;
      // else cout << "  StartVertex: ???" << endl;
      cout << "  Momentum: ("
        << recpart.getMomentum().x << ", "
        << recpart.getMomentum().y << ", "
        << recpart.getMomentum().z << ")" << endl;
      cout << "  p=|Momentum|: " << edm4hep::utils::p(recpart) << endl;
      cout << "  Energy:       " << recpart.getEnergy() << endl;
      cout << "  # of clusters: "   << recpart.clusters_size()    << endl;
      cout << "  # of tracks:   "   << recpart.tracks_size()      << endl;
      cout << "  # of PIDs:     "   << recpart.particleIDs_size() << endl;
      cout << "  # of combined reconstructed parts: " << recpart.particles_size() << endl;
      // for(auto track : recpart.getTracks()) {
      //   // ...
      // }
      // for(auto cluster : recpart.getClusters()) {
      //   // ...
      // }
    }


    // read kinematics calculations from upstream /////////////////////////
    // TODO: cross check these with our calculations from `Kinematics`
    cout << endl << "KINEMATICS, calculated from upstream:" << endl;
    std::vector<std::string> upstreamReconMethodList = {
      "Truth",
      "Electron",
      "DA",
      "JB",
      "Sigma"
    };
    printf("  %10s %8s %8s %8s %8s %8s %8s\n",
        "",
        "x",
        "Q2",
        "W",
        "y",
        "nu",
        "elec?"
        );
    for(auto upstreamReconMethod : upstreamReconMethodList)
      for(const auto& calc : evStore.get<eicd::InclusiveKinematicsCollection>("InclusiveKinematics"+upstreamReconMethod) )
        printf("  %10s %8.2f %8.2f %8.2f %8.2f %8.2f %8d\n",
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
  // event loop end =========================================================

  // finish execution
  podioReader.closeFile();
  Finish();
}
