// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Sanghwa Park, Christopher Dilks

#ifndef AnalysisAthena_
#define AnalysisAthena_

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "Analysis.h"

class Clusters
{
  public:
    Clusters() {}
    Clusters(double E_, double x_, double y_, double z_, double theta_, double phi_) {}
    virtual ~Clusters() {}

    double E;
    double x;
    double y;
    double z;
    double theta;
    double phi;

};

class Particles
{
  public:
    int pid;
    int charge;
    int mcID;
    TLorentzVector vecPart;
};

class AnalysisAthena : public Analysis
{
  public:
    AnalysisAthena(TString configFileName_="", TString outfilePrefix_="");
    ~AnalysisAthena();

    void Execute() override;

    ClassDefOverride(AnalysisAthena,1);
};

#endif
