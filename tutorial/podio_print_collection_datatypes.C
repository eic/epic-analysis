// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// print PODIO Collection data types for a ROOT file

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

void podio_print_collection_datatypes(std::string root_file_name) {

  // start ROOTFrameReader
  fmt::print("{:=^50}\nOpening file \"{}\"\n", "", root_file_name);
  auto reader = podio::ROOTFrameReader();
  reader.openFile(root_file_name);

  // list available data models
  fmt::print("{:=^50}\n", " Data Models ");
  for(auto model : reader.getAvailableDatamodels())
    fmt::print(" - {}\n", model);

  // loop over categories (trees)
  for(auto cat_name : reader.getAvailableCategories()) {
    fmt::print("{:=^50}\n", fmt::format(" Read category (tree) '{}' ", cat_name));

    // read just the first event
    auto frame      = podio::Frame(reader.readNextEntry(std::string(cat_name)));
    auto id_table   = frame.getCollectionIDTableForWrite();
    auto coll_names = id_table.names();

    // get collection types
    std::vector<std::string> print_lines;
    for(const auto coll_name : coll_names) {
      const auto coll_id    = id_table.collectionID(coll_name);
      auto coll_base_ptr    = frame.get(coll_name);
      if(coll_base_ptr!=nullptr) {
        auto type_name = coll_base_ptr->getTypeName(); // collection type
        // auto type_name = coll_base_ptr->getValueTypeName(); // object type
        // auto type_name = coll_base_ptr->getDataTypeName(); // data type
        print_lines.push_back({
            fmt::format("{:<5} {:>45} - {:<45}", coll_id, coll_name, type_name)
            });
      }
      /*
      else // benign issue, see https://github.com/eic/EICrecon/issues/643
        fmt::print(stderr,"WARNING: collection '{}' (id={}) found in collection table, but not in tree '{}'\n",
            coll_name, coll_id, cat_name);
            */
    }

    // print collection types
    fmt::print("\n{:=^100}\n", fmt::format(" Collections in category (tree) '{}' ", cat_name));
    fmt::print("{:<5} {:>45} - {:<45}\n",    "ID", "Collection Name", "Collection Type");
    fmt::print("{:-<5} {:->45} - {:-<45}\n", "",   "",                "");
    for(auto line : print_lines)
      fmt::print("{}\n", line);
    fmt::print("\n");

  }
}
