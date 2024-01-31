// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// count number of events in a given file

#include <fmt/format.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TRegexp.h>
#include <Histos.h>
#include <adage/BinSet.h>

int main(int argc, char** argv) {

  // arguments
  if(argc<=2) {
    fmt::print(stderr,"USAGE: {} [output_root_file] [input_root_files]...\n",argv[0]);
    fmt::print(stderr,"       (similar syntax as `hadd`)\n");
    return 2;
  }

  // load dictionaries
  gSystem->Load("libEpicAnalysis.so");

  // get output file and list of input files
  TString out_file_name(argv[1]);
  std::vector<TString> in_file_names;
  for(int i=2; i<argc; i++)
    in_file_names.push_back(TString(argv[i]));
  fmt::print("Running {}\n", argv[0]);
  fmt::print("\n{:=<50}\n", "Input Files ");
  for(auto in_file_name : in_file_names)
    fmt::print("{}\n", in_file_name);
  fmt::print("\n{:=<50}\n{}\n{:=<50}\n", "Output File ", out_file_name, "");

  // TFile pointers
  auto out_file = new TFile(out_file_name, "RECREATE");
  std::vector<TFile*> in_files;
  for(auto in_file_name : in_file_names)
    in_files.push_back(new TFile(in_file_name));

  // object maps
  std::map<TString,TList*> in_trees;
  std::map<TString,Histos*> in_histos;
  std::map<TString,BinSet*> in_binsets;

  // tree tracker
  std::vector<TString> tree_tracker;
    
  // Merged tree vectors
  std::vector<double> xs_total = {0};
  std::vector<double> weight_track_total = {0};
  std::vector<int>   total_events = {0};
    
  // Temporary vectors
  std::vector<double>* tmp_xs_total;
  std::vector<double>* tmp_weight_track_total;
  std::vector<int>*   tmp_total_events;
    
  // loop over input files
  bool first_file = true;
  for(auto in_file : in_files) {

    // clear ttree tracker for this file
    // helps avoid duplicate trees during merging (current/backup cycles)
    tree_tracker.clear();
      
    // loop over input file's keys
    TListIter nextKey(in_file->GetListOfKeys());
    while(TKey * key = (TKey*) nextKey()) {
      TString keyName(key->GetName());
      // handle TTrees: add their pointers to lists `in_trees`
      if(TString(key->GetClassName())=="TTree") {
        if(std::find(tree_tracker.begin(), tree_tracker.end(), keyName) != tree_tracker.end()){
            continue; // TTree named 'keyName' already added to list, this one way a backup so skip it
        }
        auto obj = (TTree*) key->ReadObj(); 
        if(first_file) {
          fmt::print("Start list for TTree '{}'\n",keyName);
          in_trees.insert({ keyName, new TList() });
        }
        try {
          in_trees.at(keyName)->Add(obj);
          tree_tracker.push_back(keyName);
        } catch(const std::out_of_range& e) {
          fmt::print(stderr,"ERROR: cannot find TTree '{}' in {}\n",keyName,in_file->GetName());
        }
      }

      // if it's a directory of Histos histograms, ignore it since we will merge the Histos objects instead
      else if(TString(key->GetClassName())=="TDirectoryFile" && keyName.Contains(TRegexp("^histArr_"))) {
        // if(first_file) fmt::print("Skip Histos TDirectory {}\n",keyName);
      }

      // handle Histos objects: add all histograms to the first file's histograms
      else if(TString(key->GetClassName())=="Histos") {
        auto obj = (Histos*) key->ReadObj(); 
        if(first_file) {
          fmt::print("Read Histos '{}'\n",keyName);
          in_histos.insert({ keyName, obj });
        } else {
          Histos *out_histos;
          try {
            out_histos = in_histos.at(keyName);
            out_histos->AddHistos(obj);
          } catch(const std::out_of_range& e) {
            fmt::print(stderr,"ERROR: cannot find Histos '{}' in {}\n",keyName,in_file->GetName());
          }
        }
      }

      // handle BinSet objects: store the first file's BinSets; TODO: compare with all the others
      else if(TString(key->GetClassName())=="BinSet") {
        if(first_file) {
          auto obj = (BinSet*) key->ReadObj(); 
          fmt::print("Read BinSet '{}'\n",keyName);
          in_binsets.insert({ keyName, obj });
        }
      }

      // handle XsTotal
      else if(keyName=="XsTotal"){
          in_file->GetObject(keyName,tmp_xs_total);
          if (tmp_xs_total->at(0)>xs_total.at(0)) // Replace XsTotal with larger one
              xs_total[0] = tmp_xs_total->at(0);
      }
        
      // handle WeightTrackTotal
      else if(keyName=="WeightTrackTotal"){
          in_file->GetObject(keyName,tmp_weight_track_total);
          weight_track_total[0]+=tmp_weight_track_total->at(0);
      }
    
      // handle TotalEvents
      else if(keyName=="TotalEvents"){
          in_file->GetObject(keyName,tmp_total_events);
          total_events[0]+=tmp_total_events->at(0);
      }
      // handle Weights // TODO
      else if(
          TString(key->GetClassName()).Contains(TRegexp("vector<.*>")) &&
          ( keyName=="XsTotal" || keyName.Contains(TRegexp("^Weight.*Total$")) )
          )
      {
        if(first_file)
          fmt::print(stderr,"WARNING: merging of {} object named '{}' is not yet supported (TODO)\n",key->GetClassName(),keyName);
      }

      // unknown object
      else {
        if(first_file)
          fmt::print(stderr,"WARNING: not sure how to handle {} object named '{}'. Ignoring.\n",key->GetClassName(),keyName);
      }

    } // end TKey loop
    first_file = false;
  } // end in_file loop

  // focus to output directory
  fmt::print("{:=<50}\n","");
  out_file->cd();

  // merge TTrees in `in_trees` and write
  for(auto [name,list] : in_trees) {
    fmt::print("Merging input trees for tree '{}'\n",name);
    TTree::MergeTrees(list)->Write(name);
  }

  // write Histos objects
  for(auto [name,histos] : in_histos) {
    fmt::print("Writing Histos '{}'\n",name);
    histos->WriteHists(out_file);
  }
  for(auto [name,histos] : in_histos)
    histos->Write();
  fmt::print(stderr,"WARNING: any Hist4D objects were not added, since Hist4D::Add must be implemented (TODO)\n");

  // write BinSets
  for(auto [name,binset] : in_binsets) {
    fmt::print("Writing BinSet '{}'\n",name);
    binset->Write(name);
  }

  // write std::vectors
  out_file->WriteObject(&xs_total,"XsTotal");
  out_file->WriteObject(&weight_track_total,"WeightTrackTotal");
  out_file->WriteObject(&total_events,"TotalEvents");
    
  // close files
  for(auto in_file : in_files)
    in_file->Close();
  out_file->Close();

  fmt::print("{:=<50}\n{} written.\n\n","",out_file_name);
  return 0;
}
