// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Sanghwa Park, Christopher Dilks

#include "AnalysisAthena.h"

using std::map;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

// constructor
AnalysisAthena::AnalysisAthena(TString configFileName_, TString outfilePrefix_) :
  Analysis(configFileName_, outfilePrefix_)
{ };

// destructor
AnalysisAthena::~AnalysisAthena() { };


//=============================================
// perform the analysis
//=============================================
void AnalysisAthena::Execute()
{
  // setup
  Prepare();

  // read dd4hep tree
  auto chain = std::make_unique<TChain>("events");
  for(Int_t idx=0; idx<infiles.size(); ++idx) {
    for(std::size_t idxF=0; idxF<infiles[idx].size(); ++idxF) {
      // std::cout << "Adding " << infiles[idx][idxF] << std::endl;
      chain->Add(infiles[idx][idxF].c_str());
    }
  }
  chain->CanDeleteRefs();

  TTreeReader tr(chain.get());

  // Truth
  TTreeReaderArray<Int_t>    mcparticles_ID(tr,        "mcparticles.ID");
  TTreeReaderArray<Int_t>    mcparticles_pdgID(tr,     "mcparticles.pdgID");
  TTreeReaderArray<Double_t> mcparticles_psx(tr,       "mcparticles.ps.x");
  TTreeReaderArray<Double_t> mcparticles_psy(tr,       "mcparticles.ps.y");
  TTreeReaderArray<Double_t> mcparticles_psz(tr,       "mcparticles.ps.z");
  TTreeReaderArray<Int_t>    mcparticles_status(tr,    "mcparticles.status");
  TTreeReaderArray<Int_t>    mcparticles_genStatus(tr, "mcparticles.genStatus");
  TTreeReaderArray<Double_t> mcparticles_mass(tr,      "mcparticles.mass");

  // Reco
  TTreeReaderArray<Int_t> ReconstructedParticles_pid(tr,   "ReconstructedParticles.pid");
  TTreeReaderArray<float> ReconstructedParticles_energy(tr,  "ReconstructedParticles.energy");
  TTreeReaderArray<float> ReconstructedParticles_p_x(tr,     "ReconstructedParticles.p.x");
  TTreeReaderArray<float> ReconstructedParticles_p_y(tr,     "ReconstructedParticles.p.y");
  TTreeReaderArray<float> ReconstructedParticles_p_z(tr,     "ReconstructedParticles.p.z");
  TTreeReaderArray<float> ReconstructedParticles_p(tr,       "ReconstructedParticles.momentum");
  TTreeReaderArray<float> ReconstructedParticles_th(tr,      "ReconstructedParticles.direction.theta");
  TTreeReaderArray<float> ReconstructedParticles_phi(tr,     "ReconstructedParticles.direction.phi");
  TTreeReaderArray<float> ReconstructedParticles_mass(tr,    "ReconstructedParticles.mass");
  TTreeReaderArray<short> ReconstructedParticles_charge(tr,  "ReconstructedParticles.charge");
  TTreeReaderArray<int>   ReconstructedParticles_mcID(tr,    "ReconstructedParticles.mcID.value");

  // calculate Q2 weights
  CalculateEventQ2Weights();

  // counters
  Long64_t numNoBeam, numEle, numNoEle, numNoHadrons, numProxMatched;
  numNoBeam = numEle = numNoEle = numNoHadrons = numProxMatched = 0;

  // event loop =========================================================
  cout << "begin event loop..." << endl;
  tr.SetEntriesRange(1,maxEvents);
  do {
    if(tr.GetCurrentEntry()%10000==0) cout << tr.GetCurrentEntry() << " events..." << endl;

    // resets
    kin->ResetHFS();
    kinTrue->ResetHFS();

    // generated truth loop
    /* - add truth particle to `mcpart`
     * - add to hadronic final state sums (momentum, sigma, etc.)
     * - find scattered electron
     * - find beam particles
     */
    std::vector<Particles> mcpart;
    double maxP = 0;
    int genEleID = -1;
    bool foundBeamElectron = false;
    bool foundBeamIon = false;
    for(int imc=0; imc<mcparticles_pdgID.GetSize(); imc++) {

      int pid_ = mcparticles_pdgID[imc];
      int genStatus_ = mcparticles_genStatus[imc]; // genStatus 4: beam particle,  1: final state
      double px_ = mcparticles_psx[imc];
      double py_ = mcparticles_psy[imc];
      double pz_ = mcparticles_psz[imc];
      double mass_ = mcparticles_mass[imc]; // in GeV
      double p_ = sqrt(pow(mcparticles_psx[imc],2) + pow(mcparticles_psy[imc],2) + pow(mcparticles_psz[imc],2));

      if(genStatus_ == 1) { // final state

        // add to `mcpart`
        Particles part;
        part.pid = pid_;
        part.vecPart.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
        part.mcID = mcparticles_ID[imc];
        mcpart.push_back(part);

        // add to hadronic final state sums
        kinTrue->AddToHFS(part.vecPart);

        // identify scattered electron by max momentum
        if(pid_ == 11) {
          if(p_ > maxP) {
            maxP = p_;
            kinTrue->vecElectron.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
            genEleID = mcparticles_ID[imc];
          }
        }
      }

      else if(genStatus_ == 4) { // beam particles
        if(pid_ == 11) { // electron beam
          if(!foundBeamElectron) {
            foundBeamElectron = true;
            kinTrue->vecEleBeam.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
          }
          else { ErrorPrint("ERROR: Found two beam electrons in one event"); }
        }
        else { // ion beam
          if(!foundBeamIon) {
            foundBeamIon = true;
            kinTrue->vecIonBeam.SetPxPyPzE(px_, py_, pz_, sqrt(p_*p_ + mass_*mass_));
          }
          else { ErrorPrint("ERROR: Found two beam ions in one event"); }
        }
      }
    } // end truth loop

    // check beam finding
    if(!foundBeamElectron || !foundBeamIon) { numNoBeam++; continue; };


    // reconstructed particles loop
    /* - add reconstructed particle to `recopart`
     * - find the scattered electron
     *
     */
    std::vector<Particles> recopart;
    int recEleFound = 0;
    for(int ireco=0; ireco<ReconstructedParticles_pid.GetSize(); ireco++) {

      int pid_ = ReconstructedParticles_pid[ireco];
      if(pid_ == 0) continue; // pid==0: reconstructed tracks with no matching truth pid

      // add reconstructed particle `part` to `recopart`
      Particles part;
      part.pid = pid_;
      part.mcID = ReconstructedParticles_mcID[ireco];
      part.charge = ReconstructedParticles_charge[ireco];
      double reco_E = ReconstructedParticles_energy[ireco];
      double reco_px = ReconstructedParticles_p_x[ireco];
      double reco_py = ReconstructedParticles_p_y[ireco];
      double reco_pz = ReconstructedParticles_p_z[ireco];
      double reco_mass = ReconstructedParticles_mass[ireco];
      double reco_p = sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz);
      part.vecPart.SetPxPyPzE(reco_px, reco_py, reco_pz, sqrt(reco_p*reco_p + reco_mass*reco_mass));

      // add to `recopart` and hadronic final state sums only if there is a matching truth particle
      if(part.mcID > 0) {
        for(auto imc : mcpart) {
          if(part.mcID == imc.mcID) {
            recopart.push_back(part);
            kin->AddToHFS(part.vecPart);
	    if( writeHFSTree ){
	      kin->AddToHFSTree(part.vecPart, part.pid);
	    }	    
            break;
          }
        }
      }

      // find scattered electron, by matching to truth // TODO: not realistic... is there an upstream electron finder?
      if(pid_ == 11 && part.mcID == genEleID) {
        recEleFound++;
        kin->vecElectron.SetPxPyPzE(reco_px, reco_py, reco_pz, sqrt(reco_p*reco_p + reco_mass*reco_mass));
      }

    } // end reco loop

    // skip the event if the scattered electron is not found and we need it
    if(recEleFound < 1) {
      numNoEle++;
      continue; // TODO: only need to skip if we are using a recon method that needs it (`if reconMethod==...`)
    }
    else if(recEleFound>1) ErrorPrint("WARNING: found more than 1 reconstructed scattered electron in an event");
    else numEle++;

    // subtract electron from hadronic final state variables
    kin->SubtractElectronFromHFS();
    kinTrue->SubtractElectronFromHFS();

    // skip the event if there are no reconstructed particles (other than the
    // electron), otherwise hadronic recon methods will fail
    if(kin->countHadrons == 0) {
      numNoHadrons++;
      continue;
    };
    
    // calculate DIS kinematics
    if(!(kin->CalculateDIS(reconMethod))) continue; // reconstructed
    if(!(kinTrue->CalculateDIS(reconMethod))) continue; // generated (truth)

    // Get the weight for this event's Q2
    auto Q2weightFactor = GetEventQ2Weight(kinTrue->Q2, inLookup[chain->GetTreeNumber()]);

    // fill inclusive histograms, if only `inclusive` is included in output
    // (otherwise they will be filled in track and jet loops)
    if(includeOutputSet["inclusive_only"]) {
      auto wInclusive = Q2weightFactor * weightInclusive->GetWeight(*kinTrue);
      wInclusiveTotal += wInclusive;
      FillHistosInclusive(wInclusive);
    }

    // loop over reconstructed particles again
    /* - calculate hadron kinematics
     * - fill output data structures (Histos, SidisTree, etc.)
     */
    for(auto part : recopart) {
      int pid_ = part.pid;
      int mcid_ = part.mcID;

      // final state cut
      // - check PID, to see if it's a final state we're interested in for
      //   histograms; if not, proceed to next track
      auto kv = PIDtoFinalState.find(pid_);
      if(kv!=PIDtoFinalState.end()) finalStateID = kv->second; else continue;
      if(activeFinalStates.find(finalStateID)==activeFinalStates.end()) continue;

      // calculate reconstructed hadron kinematics
      kin->vecHadron = part.vecPart;
      kin->CalculateHadronKinematics();

      // add selected single hadron FS to HFS tree
      if( writeHFSTree ){
        kin->AddTrackToHFSTree(part.vecPart, part.pid);
      }

      // find the matching truth hadron using mcID, and calculate its kinematics
      if(mcid_ > 0) {
        for(auto imc : mcpart) {
          if(mcid_ == imc.mcID) {
            kinTrue->vecHadron = imc.vecPart;
	    // add tracks of interest for kinematic studies to HFSTree
            if( writeHFSTree ){
              kinTrue->AddTrackToHFSTree(imc.vecPart, imc.pid);
            }
            break;
          }
        }
      }
      /* // deprecated, since existence of truth match is checked earlier; in practice prox matching was never called
      else {
        // give it another shot: proximity matching
        double mineta = 4.0;
        numProxMatched++;
        for(int imc=0; imc<(int)mcpart.size(); imc++) {
          if(pid_ == mcpart[imc].pid) {
            double deta = abs(kin->vecHadron.Eta() - mcpart[imc].vecPart.Eta());
            if(deta < mineta) {
              mineta = deta;
              kinTrue->vecHadron = mcpart[imc].vecPart;
            }
          }
        }
      }
      */
      kinTrue->CalculateHadronKinematics();

      // asymmetry injection
      // kin->InjectFakeAsymmetry(); // sets tSpin, based on reconstructed kinematics
      // kinTrue->InjectFakeAsymmetry(); // sets tSpin, based on generated kinematics
      // kin->tSpin = kinTrue->tSpin; // copy to "reconstructed" tSpin

      // weighting
      auto wTrack = Q2weightFactor * weightTrack->GetWeight(*kinTrue);
      wTrackTotal += wTrack;

      if(includeOutputSet["1h"]) {
        // fill single-hadron histograms in activated bins
        FillHistos1h(wTrack);
        FillHistosInclusive(wTrack);

        // fill simple tree
        // - not binned
        // - `IsActiveEvent()` is only true if at least one bin gets filled for this track
        if( writeSidisTree && HD->IsActiveEvent() ) ST->FillTree(wTrack);
      }

    }//hadron loop
    
    // fill HFSTree for each event
    if( writeHFSTree && kin->nHFS > 0) HFST->FillTree(Q2weightFactor);
    
  } while(tr.Next()); // tree reader loop
  cout << "end event loop" << endl;
  // event loop end =========================================================

  // finish execution
  Finish();

  // final printout
  cout << "Total number of scattered electrons found: " << numEle << endl;
  if(numNoEle>0)
    cerr << "WARNING: skipped " << numNoEle << " events which had no reconstructed scattered electron" << endl;
  if(numNoHadrons>0)
    cerr << "WARNING: skipped " << numNoHadrons << " events which had no reconstructed hadrons" << endl;
  if(numNoBeam>0)
    cerr << "WARNING: skipped " << numNoBeam << " events which had no beam particles" << endl;
  if(numProxMatched>0)
    cerr << "WARNING: " << numProxMatched << " recon. particles were proximity matched to truth (when mcID match failed)" << endl;

}
