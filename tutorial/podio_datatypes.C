// print PODIO data types for a ROOT file
//
// USAGE:
// root -b -q tutorial/podio_datatypes.C'("","my_root_file.root")'  // print all collections
// root -b -q tutorial/podio_datatypes.C'("MCParticles","my_root_file.root")'  // print MCParticles only
//

R__LOAD_LIBRARY(podioDict)
R__LOAD_LIBRARY(podioRootIO)
R__LOAD_LIBRARY(edm4hep)
R__LOAD_LIBRARY(edm4eic)
// R__LOAD_LIBRARY(eicd) // replaced with edm4eic
R__LOAD_LIBRARY(fmt)
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"
#include "podio/CollectionBase.h"
#include "fmt/format.h"

void podio_datatypes(
    std::string branch_name    = "", // leave empty to print all
    const char *root_file_name = "datarec/epic_test/rec-dis_10x100_minQ2=1.root"
    )
{
  // open event store
  podio::ROOTReader reader;
  podio::EventStore store;
  reader.openFile(root_file_name);
  store.setReader(&reader);

  // print header
  fmt::print("\n{:+^55}\n"," podio collection id table ");
  fmt::print("\n{:<4} {:<50}\n", "ID", "CollectionName");
  fmt::print("{:>50}\n{:>50}\n{:>50}\n",
      "CollectionTypeName", 
      "ValueTypeName",
      "DataTypeName"
      );

  // get list of collection names
  const auto collIDTable = store.getCollectionIDTable();
  auto collNames = collIDTable->names();
  if(branch_name!="") {
    collNames.clear();
    collNames.push_back(branch_name);
  }
  // std::vector<std::string> collNames = {
  //   "MCParticles",
  //   "ReconstructedParticles",
  //   "ReconstructedParticlesAssoc"
  // };

  // print collection types
  for(const auto collName : collNames) {
    const auto collID = collIDTable->collectionID(collName);
    fmt::print("\n{:<4} {:<50}\n",collID,collName);
    podio::CollectionBase *collBase;
    if(store.get(collID, collBase))
      fmt::print("{:>50}\n{:>50}\n{:>50}\n",
          collBase->getTypeName(),
          collBase->getValueTypeName(),
          collBase->getDataTypeName()
          );
    else
      fmt::print("{:>50}\n","not found in event store");
  }
  fmt::print("\n{:+^55}\n","");
}
