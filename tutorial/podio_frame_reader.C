// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// example how to read a ROOT file with `PODIO` `ROOTFrameReader`

R__LOAD_LIBRARY(fmt)
#include "fmt/format.h"

R__LOAD_LIBRARY(podio)
R__LOAD_LIBRARY(podioRootIO)
#include "podio/ROOTFrameReader.h"
#include "podio/Frame.h"

R__LOAD_LIBRARY(edm4hep)
R__LOAD_LIBRARY(edm4eic)
R__LOAD_LIBRARY(edm4hepDict)
R__LOAD_LIBRARY(edm4eicDict)
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4hep/utils/kinematics.h"

// -----------------------------------------------------------------------------

void PrintParticle(const edm4hep::MCParticle& P);
void PrintParticle(const edm4eic::ReconstructedParticle& P);

// -----------------------------------------------------------------------------

void podio_frame_reader(
    const std::string& root_file_name, // input ROOT file
    unsigned max_num_entries = 0       // maximum number of entries to read
    )
{

  // start ROOTFrameReader
  fmt::print("{:=^50}\nOpening file \"{}\"\n", "", root_file_name);
  auto reader = podio::ROOTFrameReader();
  reader.openFile(root_file_name);

  // get number of events
  const std::string tree_name = "events"; // could be obtained from `reader.getAvailableCategories()`
  auto num_entries = reader.getEntries(tree_name);
  fmt::print("Number of entries in \"{}\": {}\n", tree_name, num_entries);
  if(max_num_entries>0) num_entries = std::min(num_entries, max_num_entries);
  fmt::print("Reading {} of them...\n", num_entries);

  // event loop
  for(unsigned e=0; e<num_entries; e++) {

    // next event
    auto frame = podio::Frame(reader.readNextEntry(tree_name));

    // dump list of collections (from first event)
    if(e==0) {
      fmt::print("\n{:=^50}\n", " COLLECTIONS ");
      for(auto coll : frame.getAvailableCollections())
        fmt::print("  {}\n", coll);
      fmt::print("\n{:=^50}\n", " BEGIN EVENT LOOP ");
    }

    // print event number
    fmt::print("\n{:=<50}\n", fmt::format("=== EVENT {} ",e));

    // get collections
    std::string rec_parts_collname  = "ReconstructedChargedParticles";
    std::string rec_assocs_collname = "ReconstructedChargedParticlesAssociations";
    auto& rec_parts  = frame.get<edm4eic::ReconstructedParticleCollection>(rec_parts_collname);
    auto& rec_assocs = frame.get<edm4eic::MCRecoParticleAssociationCollection>(rec_assocs_collname);

    // check for existence of these collections
    if( ! (rec_parts.isValid() && rec_assocs.isValid())) {
      fmt::print(stderr,"WARNING: missing collections for this event\n");
      continue;
    }

    // get the number of entries in each collection
    fmt::print("Number of entries in collections:\n");
    fmt::print(" {:>50}: {}\n", rec_parts_collname,  rec_parts.size());
    fmt::print(" {:>50}: {}\n", rec_assocs_collname, rec_assocs.size());

    // warn if number of reconstructed particles != number of associations
    // - you may not want this warning or check for this, because sometimes not
    //   every reconstructed particle will have an associated MC particle
    if(rec_parts.size() != rec_assocs.size())
      fmt::print(stderr,"WARNING: {} rec. particles != {} rec. particle associations\n", rec_parts.size(), rec_assocs.size());

    // loop over reconstructed particles
    /*
    fmt::print("\n{:-^50}\n", fmt::format("Reconstructed Particles ",e));
    for(const auto& rec_part : rec_parts) {
      PrintParticle(rec_part);
      fmt::print("\n");
    }
    */

    // loop over reconstructed particle associations
    fmt::print("\n{:-^50}\n", fmt::format("Reconstructed Particle Associations ",e));
    for(const auto& rec_assoc : rec_assocs) {
      PrintParticle(rec_assoc.getRec());
      PrintParticle(rec_assoc.getSim());
      fmt::print("\n");
    }

  } // end event loop

}

// -----------------------------------------------------------------------------

void PrintParticle(const edm4hep::MCParticle& P) { 
  fmt::print("=> True MC Particle\n");
  fmt::print("  {:>20}: {}\n",     "True PDG",     P.getPDG()             );
  fmt::print("  {:>20}: {}\n",     "Status",       P.getGeneratorStatus() );
  fmt::print("  {:>20}: {} GeV\n", "Energy",       P.getEnergy()          );
  fmt::print("  {:>20}: {} GeV\n", "p=|Momentum|", edm4hep::utils::p(P)   );
  fmt::print("  {:>20}: {} GeV\n", "pT_lab",       edm4hep::utils::pT(P)  );
  fmt::print("  {:>20}: ({}, {}, {}) GeV\n",
      "3-Momentum",
      P.getMomentum().x,
      P.getMomentum().y,
      P.getMomentum().z
      );
  fmt::print("  {:>20}: ({}, {}, {}) mm\n",
      "Vertex",
      P.getVertex().x,
      P.getVertex().y,
      P.getVertex().z
      );
  fmt::print("  {:>20}:\n", "Parents");
  if(P.parents_size()>0) {
    for(const auto& parent : P.getParents())
      fmt::print("    {:>20}: {}\n", "PDG", parent.getPDG());
  } else fmt::print("    {:>20}\n", "None");
  fmt::print("  {:>20}:\n", "Daughters");
  if(P.daughters_size()>0) {
    for(const auto& daughter : P.getDaughters())
      fmt::print("    {:>20}: {}\n", "PDG", daughter.getPDG());
  } else fmt::print("    {:>20}\n", "None");
}

// -----------------------------------------------------------------------------

void PrintParticle(const edm4eic::ReconstructedParticle& P) {
  fmt::print("=> Reconstructed Particle:\n");
  fmt::print("  {:>20}: ", "ParticleIDUsed::PDG");
  if(P.getParticleIDUsed().isAvailable()) fmt::print("{}\n", P.getParticleIDUsed().getPDG());
  else fmt::print("???\n");
  fmt::print("  {:>20}: {} GeV\n", "Mass",         P.getMass()           );
  fmt::print("  {:>20}: {}\n",     "Charge",       P.getCharge()         );
  fmt::print("  {:>20}: {} GeV\n", "Energy",       P.getEnergy()         );
  fmt::print("  {:>20}: {} GeV\n", "p=|Momentum|", edm4hep::utils::p(P)  );
  fmt::print("  {:>20}: {} GeV\n", "pT_lab",       edm4hep::utils::pT(P) );
  fmt::print("  {:>20}: ({}, {}, {}) GeV\n",
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
