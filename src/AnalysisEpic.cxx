#include "AnalysisEpic.h"

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
  fmt::print("begin event loop...\n");
  for(unsigned e=0; e<ENT; e++) {
    if(e%10000==0) fmt::print("{} events...\n",e);
    if(verbose) fmt::print("\n\n{:=<70}\n",fmt::format("EVENT {} ",e));

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();
    double mcPartElectronP   = 0.0;
    bool double_counted_beam = false;
    int num_ele_beams        = 0;
    int num_ion_beams        = 0;
    int num_sim_electrons    = 0;
    int num_rec_electrons    = 0;
    
    // read particle collections for this event
    // FIXME: not yet fully using `edm4*` in physics_benchmarks pipelines; instead using `eicd`
    const auto& simParts   = evStore.get<edm4hep::MCParticleCollection>("MCParticles");
    const auto& recParts   = evStore.get<eicd::ReconstructedParticleCollection>("ReconstructedParticles");
    const auto& mcRecLinks = evStore.get<eicd::MCRecoParticleAssociationCollection>("ReconstructedParticlesAssoc");

    // data objects
    edm4hep::MCParticle mcPartEleBeam;
    edm4hep::MCParticle mcPartIonBeam;
    edm4hep::MCParticle mcPartElectron;
    std::set<eicd::ReconstructedParticle> recPartsToAnalyze;

    // loop over generated particles
    if(verbose) fmt::print("\n{:-<60}\n","MCParticles ");
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
            ErrorPrint(fmt::format("WARNING: Unknown beam particle with PDG={}",pid));
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

    } // end loop over generated particles

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
      if(verbose) fmt::print("\n{:-<60}\n","GENERATED BEAMS ");
      PrintParticle(mcPartEleBeam);
      PrintParticle(mcPartIonBeam);
      if(verbose) fmt::print("\n{:-<60}\n","GENERATED SCATTERED ELECTRON ");
      PrintParticle(mcPartElectron);
    }

    // loop over reconstructed particles
    /*
    if(verbose) fmt::print("\n{:-<60}\n","ReconstructedParticles ");
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
    if(verbose) fmt::print("\n{:-<60}\n","MC<->Reco ASSOCIATIONS ");
    for(const auto& link : mcRecLinks ) {
      auto recPart = link.getRec();
      auto simPart = link.getSim();
      bool truthMatch = simPart.isAvailable();
      // if(!truthMatch) continue; // FIXME: consider using this once we have matching

      // print out this reconstructed particle, and its match
      if(verbose) {
        fmt::print("\n   {:->35}\n"," reconstructed particle:");
        PrintParticle(recPart);
        fmt::print("\n   {:.>35}\n"," truth match:");
        if(truthMatch) PrintParticle(simPart);
        else fmt::print("     {:>35}\n","NO MATCH");
        fmt::print("\n");
      }

      // get PID
      bool usedTruthPID = false;
      auto pid = GetReconstructedPDG(simPart, recPart, usedTruthPID);
      if(verbose) fmt::print("   GetReconstructedPDG = {}\n",pid);
      // if(usedTruthPID) continue; // FIXME: consider using this once we have decent PID

      // add to list of reconstructed particles to analyze, and to the HFS
      recPartsToAnalyze.insert(recPart);
      kin->AddToHFS(GetP4(recPart));

      // find scattered electron, by matching to truth
      // FIXME: not working unless we have truth matching and/or reconstructed PID
      // FIXME: does `simPart==mcPartElectron` actually work !? - alternatively use ID to check matching
      /*
      if(pid==constants::pdgElectron && simPart==mcPartElectron) {
        num_rec_electrons++;
        kin->vecElectron = GetP4(recPart);
      }
      */

    } // end loop over MC<->Rec associations

    /* // FIXME: beyond here, need scattered electron
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

    //// TODO: stopped syncing with AnalysisAthena here ////

    */





    // read kinematics calculations from upstream /////////////////////////
    // TODO: cross check these with our calculations from `Kinematics`
    fmt::print("\n{:-<60}\n","KINEMATICS, calculated from upstream: ");
    fmt::print("  {:>10} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\n", "", "x", "Q2", "W", "y", "nu", "elec?");
    for(const auto upstreamReconMethod : upstreamReconMethodList)
      for(const auto& calc : evStore.get<eicd::InclusiveKinematicsCollection>("InclusiveKinematics"+upstreamReconMethod) )
        fmt::print("  {:10} {:8.5f} {:8.2f} {:8.2f} {:8.5f} {:8.2f} {:>8}\n",
            upstreamReconMethod,
            calc.getX(),
            calc.getQ2(),
            calc.getW(),
            calc.getY(),
            calc.getNu(),
            calc.getScat().isAvailable()
            );


    // next event //////////////////////////////////
    evStore.clear();
    podioReader.endOfEvent();

  } // event loop
  fmt::print("end event loop\n");

  // finish execution
  podioReader.closeFile();
  Finish();
}


// particle printers //////////////////////////////////////////////

void AnalysisEpic::PrintParticle(const edm4hep::MCParticle& P) { 
  fmt::print("\n");
  fmt::print("  {:>20}: {}\n", "PDG",          P.getPDG()             );
  fmt::print("  {:>20}: {}\n", "Status",       P.getGeneratorStatus() );
  fmt::print("  {:>20}: {}\n", "Energy",       P.getEnergy()          );
  fmt::print("  {:>20}: {}\n", "p=|Momentum|", edm4hep::utils::p(P)   );
  fmt::print("  {:>20}: {}\n", "pT_lab",       edm4hep::utils::pT(P)  );
  fmt::print("  {:>20}: ({}, {}, {})\n",
      "3-Momentum",
      P.getMomentum().x,
      P.getMomentum().y,
      P.getMomentum().z
      );
  fmt::print("  {:>20}: ({}, {}, {})\n",
      "Vertex",
      P.getVertex().x,
      P.getVertex().y,
      P.getVertex().z
      );
  fmt::print("  {:>20}:\n", "Parents");
  for(const auto& parent : P.getParents())
    fmt::print("    {:>20}: {}\n", "PDG", parent.getPDG());
  fmt::print("  {:>20}:\n", "Daughters");
  for(const auto& daughter : P.getDaughters())
    fmt::print("    {:>20}: {}\n", "PDG", daughter.getPDG());
}

void AnalysisEpic::PrintParticle(const eicd::ReconstructedParticle& P) {
  fmt::print("\n");
  fmt::print("  {:>20}: ", "PDG");
  if(P.getParticleIDUsed().isAvailable()) fmt::print("{}\n", P.getParticleIDUsed().getPDG());
  else fmt::print("???\n");
  fmt::print("  {:>20}: {}\n", "Mass",         P.getMass()           );
  fmt::print("  {:>20}: {}\n", "Charge",       P.getCharge()         );
  fmt::print("  {:>20}: {}\n", "Energy",       P.getEnergy()         );
  fmt::print("  {:>20}: {}\n", "p=|Momentum|", edm4hep::utils::p(P)  );
  fmt::print("  {:>20}: {}\n", "pT_lab",       edm4hep::utils::pT(P) );
  fmt::print("  {:>20}: ({}, {}, {})\n",
      "3-Momentum",
      P.getMomentum().x,
      P.getMomentum().y,
      P.getMomentum().z
      );
  fmt::print("  {:>20}: {}\n", "# of clusters", P.clusters_size()    );
  fmt::print("  {:>20}: {}\n", "# of tracks",   P.tracks_size()      );
  fmt::print("  {:>20}: {}\n", "# of PIDs",     P.particleIDs_size() );
  fmt::print("  {:>20}: {}\n", "# of recParts", P.particles_size()   );
  // for(const auto& track : P.getTracks()) {
  //   // ...
  // }
  // for(const auto& cluster : P.getClusters()) {
  //   // ...
  // }
}


// helper methods //////////////////////////////////////////////

// get PDG from reconstructed particle
int AnalysisEpic::GetReconstructedPDG(
    const edm4hep::MCParticle& simPart,
    const eicd::ReconstructedParticle& recPart,
    bool& usedTruth
    )
{
  int pid = 0;
  usedTruth = false;

  // if using edm4hep::ReconstructedParticle:
  /*
  if(recPart.getParticleIDUsed().isAvailable()) // FIXME: not available
    pid = recPart.getParticleIDUsed().getPDG();
  */
  
  // if using eicd::ReconstructedParticle:
  // pid = recPart.getPDG(); // FIXME: not available either

  // if reconstructed PID is unavailable, use MC PDG
  if(pid==0) {
    usedTruth = true;
    if(simPart.isAvailable())
      pid = simPart.getPDG();
  }

  return pid;
}
