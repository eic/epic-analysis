// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

#include "AnalysisEpicPodio.h"

AnalysisEpicPodio::AnalysisEpicPodio(TString infileName_, TString outfilePrefix_)
  : Analysis(infileName_, outfilePrefix_)
  , crossCheckKinematics(false)
{}

AnalysisEpicPodio::~AnalysisEpicPodio() {}

void AnalysisEpicPodio::Execute()
{
  // setup
  Prepare();

  // produce flat list of files from `infiles`
  std::vector<std::string> infilesFlat;
  for(const auto fileList : infiles)
    for(const auto fileName : fileList)
      infilesFlat.push_back(fileName);

  // create PODIO reader
  podioReader.openFiles(infilesFlat);
  ENT = podioReader.getEntries(inputTreeName);
  if(maxEvents>0) ENT = std::min(maxEvents,ENT);
  
  // calculate Q2 weights
  CalculateEventQ2Weights();

  // upstream reconstruction methods
  // - list of upstream methods
  const std::vector<std::string> upstreamReconMethodList = {
    "Truth",
    "Electron",
    "DA",
    "JB",
    "Sigma"
  };
  // - association of `Kinematics::CalculateDIS` reconstruction method with upstream;
  //   for those unavailable upstream, use `"NONE"`
  const std::map<TString,std::string> associatedUpstreamMethodMap = {
    { "Ele",    "Electron" },
    { "DA",     "DA"       },
    { "JB",     "JB"       },
    { "Sigma",  "Sigma"    },
    { "Mixed",  "NONE"     },
    { "eSigma", "NONE"     }
  };
  // - get upstream method associated with `reconMethod`
  const auto& associatedUpstreamMethod = associatedUpstreamMethodMap.at(reconMethod);

  // event loop =========================================================
  fmt::print("begin event loop...\n");
  for(unsigned e=0; e<ENT; e++) {
    if(e%10000==0) fmt::print("{} events...\n",e);
    if(verbose) fmt::print("\n\n{:=<70}\n",fmt::format("EVENT {} ",e));

    // next event
    auto podioFrame = podio::Frame(podioReader.readNextEntry(inputTreeName));

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
    // FIXME: which collections to read: "ReconstructedParticles*" or "ReconstructedChargedParticles*"?
    const auto& simParts    = podioFrame.get<edm4hep::MCParticleCollection>("MCParticles");
    const auto& recParts    = podioFrame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
    const auto& mcRecAssocs = podioFrame.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedChargedParticlesAssociations");

    // data objects
    edm4hep::MCParticle mcPartEleBeam;
    edm4hep::MCParticle mcPartIonBeam;
    edm4hep::MCParticle mcPartElectron;

    // loop over generated particles
    if(verbose) fmt::print("\n{:-<60}\n","MCParticles ");
    for(const auto& simPart : simParts) {

      // print out this MCParticle
      // if(verbose) PrintParticle(simPart);

      // generated particle properties
      auto simPDG = simPart.getPDG();

      // add to Hadronic Final State (HFS) sums
      kinTrue->AddToHFS(GetP4(simPart));

      // filter for beam particles
      if(simPart.getGeneratorStatus() == constants::statusBeam) {
        switch(simPDG) {
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
            ErrorPrint(fmt::format("WARNING: Unknown beam particle with PDG={}",simPDG));
        }
      }

      // filter for scattered electron: select the one with the highest |p|
      if(simPart.getGeneratorStatus() == constants::statusFinal) {
        if(simPDG == constants::pdgElectron) {
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
    if(num_ele_beams==0)     { ErrorPrint("WARNING: missing MC electron beam");      continue; }
    if(num_ion_beams==0)     { ErrorPrint("WARNING: missing MC ion beam");           continue; }
    if(num_sim_electrons==0) { ErrorPrint("WARNING: missing scattered electron");    continue; }
    if(double_counted_beam)  { ErrorPrint("WARNING: found multiple beam particles"); continue; }

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

    // add reconstructed particles to Hadronic Final State (HFS)
    /* the following will run loops over Reconstructed Particle <-> MC Particle associations
     * - uses high-order function `LoopMCRecoAssocs` for common tasks, such as quality cuts
     *   and getting the reconstructed PID (PDG)
     * - first define a first-order function (`payload`), then call `LoopMCRecoAssocs`
     * - see `LoopMCRecoAssocs` for `payload` signature
     */
    auto AllRecPartsToHFS = [&] (const auto& simPart, const auto& recPart, auto recPDG) {
      kin->AddToHFS(GetP4(recPart));
    };
    if(verbose) fmt::print("\n{:-<60}\n","MC<->Reco ASSOCIATIONS ");
    LoopMCRecoAssocs(recParts, mcRecAssocs, AllRecPartsToHFS, verbose);

    // find reconstructed electron
    // ============================================================================
    /* FIXME: need realistic electron finder; all of the following options rely
     * on MC-truth matching; is there any common upstream realistic electron finder?
     */

    // find scattered electron by simply matching to truth
    auto FindRecoEleByTruth = [&] (const auto& simPart, const auto& recPart, auto recPDG) {
      if(recPDG==constants::pdgElectron && simPart.id()==mcPartElectron.id()) {
        num_rec_electrons++;
        kin->vecElectron = GetP4(recPart);
      }
    };
    LoopMCRecoAssocs(recParts, mcRecAssocs, FindRecoEleByTruth);

    // Fallback: use electron finder from upstream algorithm `InclusiveKinematics*`
    if(num_rec_electrons == 0) {
      ErrorPrint("WARNING: first attempt to find reconstructed electron failed; trying InclusiveKinematics collection");
      const auto& disCalcs = podioFrame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsElectron");
      if(disCalcs.size() != 1) ErrorPrint(fmt::format("WARNING: disCalcs.size = {} != 1 for this event",disCalcs.size()));
      for(const auto& calc : disCalcs) {
        auto ele = calc.getScat();
        if( ! ele.isAvailable()) {
          ErrorPrint("WARNING: `disCalcs` scattered electron unavailable");
          continue;
        }
        num_rec_electrons++;
        kin->vecElectron = GetP4(ele);
      }
    }

    // check for found reconstructed scattered electron
    if(num_rec_electrons == 0) { ErrorPrint("WARNING: reconstructed scattered electron not found"); continue; }
    if(num_rec_electrons >  1) { ErrorPrint("WARNING: found more than 1 reconstructed scattered electron"); }

    // subtract electron from hadronic final state variables
    kin->SubtractElectronFromHFS();
    kinTrue->SubtractElectronFromHFS();

    // skip the event if there are no reconstructed particles (other than the
    // electron), otherwise hadronic recon methods will fail
    if(kin->countHadrons == 0) { ErrorPrint("WARNING: no hadrons"); }
  
    // calculate DIS kinematics; skip the event if the calculation did not go well
    if( ! kin->CalculateDIS(reconMethod)     ) continue; // reconstructed
    if( ! kinTrue->CalculateDIS(reconMethod) ) continue; // generated (truth)

    // Get the weight for this event's Q2
    //   FIXME: we are in a PODIO reader event loop, thus we need an
    //          alternative to `chain->GetTreeNumber()`; currently disabling weighting
    //          for now, by setting `wTrack=1.0`
    // auto Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);
    auto Q2weightFactor = 1.0;

    // fill inclusive histograms, if only `inclusive` is included in output
    // (otherwise they will be filled in track and jet loops)
    if(includeOutputSet["inclusive_only"]) {
      auto wInclusive = Q2weightFactor * weightInclusive->GetWeight(*kinTrue);
      wInclusiveTotal += wInclusive;
      FillHistosInclusive(wInclusive);
    }

    // loop over Reco<->MC associations again
    /* - calculate SIDIS kinematics
     * - fill output data structures
     */
    auto SidisOutput = [&] (const auto& simPart, const auto& recPart, auto recPDG) {

      // final state cut
      // - check PID, to see if it's a final state we're interested in
      auto kv = PIDtoFinalState.find(recPDG);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else return;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) return;

      // set SIDIS particle 4-momenta, and calculate their kinematics
      kinTrue->vecHadron = GetP4(simPart);
      kinTrue->CalculateHadronKinematics();
      kin->vecHadron = GetP4(recPart);
      kin->CalculateHadronKinematics();

      // weighting
      auto wTrack = Q2weightFactor * weightTrack->GetWeight(*kinTrue);
      wTrackTotal += wTrack;

      // fill single-hadron histograms in activated bins
      FillHistos1h(wTrack);
      FillHistosInclusive(wTrack);

      // fill simple tree
      // - not binned
      // - `IsActiveEvent()` is only true if at least one bin gets filled for this track
      if( writeSidisTree && HD->IsActiveEvent() ) ST->FillTree(wTrack);
    };
    LoopMCRecoAssocs(recParts, mcRecAssocs, SidisOutput);


    // read kinematics calculations from upstream /////////////////////////
    // TODO: cross check these with our calculations from `Kinematics`
    if(crossCheckKinematics) {
      auto PrintRow = [] <typename T> (std::string name, std::vector<T> vals, bool header=false) {
        fmt::print("  {:>16}",name);
        if(header) { for(auto val : vals) fmt::print(" {:>8}",   val); }
        else       { for(auto val : vals) fmt::print(" {:8.4f}", val); }
        fmt::print("\n");
      };
      // upstream calculations
      fmt::print("\n{:-<75}\n","KINEMATICS, calculated from upstream: ");
      PrintRow("", std::vector<std::string>({ "x", "Q2", "W", "y", "nu" }), true);
      for(const auto upstreamReconMethod : upstreamReconMethodList)
        for(const auto& calc : podioFrame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematics"+upstreamReconMethod) )
          PrintRow( upstreamReconMethod, std::vector<float>({
              calc.getX(),
              calc.getQ2(),
              calc.getW(),
              calc.getY(),
              calc.getNu()
              }));
      // local calculations
      fmt::print("{:-<75}\n",fmt::format("KINEMATICS, calculated locally in EPIC-ANALYSIS, with method \"{}\": ",reconMethod));
      auto PrintKinematics = [&PrintRow] (std::string name, std::shared_ptr<Kinematics> K) {
        PrintRow( name, std::vector<Double_t>({
            K->x,
            K->Q2,
            K->W,
            K->y,
            K->Nu
            }));
      };
      PrintKinematics("Truth",kinTrue);
      PrintKinematics("Reconstructed",kin);
      // compare upstream and local
      if(associatedUpstreamMethod != "NONE") {
        fmt::print("{:-<75}\n",fmt::format("DIFFERENCE: upstream({}) - local({}): ",associatedUpstreamMethod,reconMethod));
        for(const auto upstreamMethod : std::vector<std::string>({"Truth",associatedUpstreamMethod})) {
          const auto& upstreamCalcs = podioFrame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematics"+upstreamMethod);
          for(const auto& upstreamCalc : upstreamCalcs) {
            auto K    = upstreamMethod=="Truth" ? kinTrue : kin;
            auto name = upstreamMethod=="Truth" ? "Truth" : "Reconstructed";
            PrintRow( name, std::vector<Double_t>({
                upstreamCalc.getX()  - K->x,
                upstreamCalc.getQ2() - K->Q2,
                upstreamCalc.getW()  - K->W,
                upstreamCalc.getY()  - K->y,
                upstreamCalc.getNu() - K->Nu
                }));
          }
        }
      }
      else fmt::print("{:-<75}\n  method \"{}\" is not available upstream\n","DIFFERENCE: ",reconMethod);
    } // if crossCheckKinematics

  } // event loop
  fmt::print("end event loop\n");

  // finish execution
  Finish();
}


// particle printers //////////////////////////////////////////////

void AnalysisEpicPodio::PrintParticle(const edm4hep::MCParticle& P) { 
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

void AnalysisEpicPodio::PrintParticle(const edm4eic::ReconstructedParticle& P) {
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


// helper methods /////////////////////////////////////////////////////////////

// common loop over Reconstructed Particle <-> MC Particle associations
/* - get PID
 * - basic quality cuts
 * - execute `payload`
//     payload signature: (simPart, recPart, reconstructed PDG)
 */
void AnalysisEpicPodio::LoopMCRecoAssocs(
    const edm4eic::ReconstructedParticleCollection&     recParts,
    const edm4eic::MCRecoParticleAssociationCollection& mcRecAssocs,
    std::function<void(const edm4hep::MCParticle&, const edm4eic::ReconstructedParticle&, int)> payload,
    bool printParticles
    )
{

  // check collection sizes equality
  /* 
   * NOTE: since it's possible that `recParts.size() != `mcRecAssocs.size()`,
   *       we first loop over `recParts`, and find the association in `mcRecAssocs`
   *       which contains it; this allows us to handle the case when we have a
   *       reconstructed particle with no associated MC truth particle
   *
   * FIXME: for now we just warn about these unequal-size cases; we should think of a way to
   *        properly handle them
   */
  if(recParts.size() > mcRecAssocs.size())
    ErrorPrint("WARNING: there are more reconstructed particles than MC-Reco associations");
  else if(recParts.size() < mcRecAssocs.size())
    ErrorPrint("ERROR: there are more MC-Reco associations than reconstructed particles; this could be an upstream issue");

  // loop over reconstructed particles
  for(const auto& recPart : recParts) {

    // find the associated MC particle
    bool assocFound = false;
    for(const auto& assoc : mcRecAssocs ) {
      if(recPart.id() == assoc.getRec().id()) {
        assocFound = true;

        // get MC truth particle
        const auto& simPart = assoc.getSim();

        // print out this reconstructed particle, and its matching truth 
        if(printParticles) {
          fmt::print("\n   {:->35}\n"," reconstructed particle:");
          PrintParticle(recPart);
          fmt::print("\n   {:.>35}\n"," truth match:");
          if(simPart.isAvailable())
            PrintParticle(simPart);
          else
            fmt::print("     {:>35}\n","NO MATCH");
          fmt::print("\n");
        }

        // handle case of missing MC truth particle
        // FIXME: handle case this better, instead of just ignoring it
        if(!simPart.isAvailable()) {
          ErrorPrint("WARNING: found reconstructed particle with MC-Reco association, but no simulated particle is linked");
          break;
        }

        // get reconstructed PDG from PID
        bool usedTruthPID = false;
        auto recPDG = GetPDG(simPart, recPart, usedTruthPID);
        if(verbose) fmt::print("   GetPDG = {}\n",recPDG);
        // if(usedTruthPID) continue; // FIXME: consider using this once we have decent PID

        // run payload
        payload(simPart, recPart, recPDG);

      }
    } // end loop over `mcRecAssocs`

    // handle reconstructed particles with no MC association
    // FIXME: what should we do about these?
    if(!assocFound) {
      ErrorPrint("WARNING: found reconstructed particle with no associated MC particle");
    }

  } // end loop over `recParts`

} // end LoopMCRecoAssocs


// get PDG from reconstructed particle; resort to true PDG, if
// PID is unavailable (sets `usedTruth` to true)
int AnalysisEpicPodio::GetPDG(
    const edm4hep::MCParticle& simPart,
    const edm4eic::ReconstructedParticle& recPart,
    bool& usedTruth
    )
{

  // get PID PDG
  usedTruth        = false;
  auto pidHypBest  = recPart.getParticleIDUsed(); // presumably the best PID hypothesis
  auto pidGoodness = recPart.getGoodnessOfPID(); // FIXME: should we use this?
  if(pidHypBest.isAvailable())
    return pidHypBest.getPDG();

  // FIXME: should we check other PID hypotheses?
  // for(auto pidHyp : recPart.getParticleIDs()) {
  //   auto pdg        = pidHyp.getPDG();
  //   auto likelihood = pidHyp.getLikelihood();
  //   // ...
  // }

  // otherwise get true PDG from associated MC particle
  usedTruth = true;
  if(simPart.isAvailable())
    return simPart.getPDG();

  // last chance // FIXME: someday we may remove `ReconstructedParticle::PDG` from EDM4eic
  ErrorPrint("WARNING: falling back to `ReconstructedParticle::PDG` member, which may be deprecated soon");
  return recPart.getPDG();
}
