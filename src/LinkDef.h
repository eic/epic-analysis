// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks, Duane Byer

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// histograms
#pragma link C++ class Hist4D+;
#pragma link C++ class Histos+;
#pragma link C++ class HistConfig+;
#pragma link C++ class Adage<Histos>+;
#pragma link C++ class HistosDAG+;

// analysis objects
#pragma link C++ class Kinematics+;
#ifndef EXCLUDE_DELPHES
#pragma link C++ class KinematicsJets+;
#endif
#pragma link C++ class SidisTree+;
#pragma link C++ class HFSTree+;
#pragma link C++ class ParticleTree+;
#pragma link C++ class Weights+;
#pragma link C++ class WeightsUniform+;
#pragma link C++ class WeightsSivers+;
#pragma link C++ class WeightsCollins+;
#pragma link C++ class WeightsProduct+;
#pragma link C++ class WeightsSum+;

// analysis event loop classes
#pragma link C++ class Analysis+;
#pragma link C++ class AnalysisEpic+;
#pragma link C++ class AnalysisAthena+;
#pragma link C++ class AnalysisEcce+;

#ifdef INCLUDE_PODIO
#pragma link C++ class AnalysisEpicPodio+;
#endif

#ifndef EXCLUDE_DELPHES
#pragma link C++ class AnalysisDelphes+;
#endif

// postprocessing
#pragma link C++ class PostProcessor+;

#endif
