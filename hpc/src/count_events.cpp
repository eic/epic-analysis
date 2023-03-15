// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// count number of events in a given file

#include <fmt/format.h>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char** argv) {
  if(argc<=1) {
    fmt::print(stderr,"USAGE: {} [root_file]\n",argv[0]);
    return 2;
  }

  auto root_file_name = std::string(argv[1]);
  auto root_file = TFile::Open(root_file_name.c_str());

  if (root_file==nullptr || root_file->IsZombie()) {
    fmt::print(stderr,"ERROR: Couldn't open input file '{}'\n",root_file_name);
    return 1;
  }

  auto tree = root_file->Get<TTree>("Delphes");                   // fastsim
  if(tree == nullptr) tree = root_file->Get<TTree>("events");     // ePIC, ATHENA
  if(tree == nullptr) tree = root_file->Get<TTree>("event_tree"); // ECCE
  if(tree == nullptr) {
    fmt::print(stderr,"ERROR: Couldn't find tree in file '{}'\n",root_file_name);
    return 1;
  }

  fmt::print("{}\n", tree->GetEntries());
  root_file->Close();
  delete root_file;
  return 0;
}
