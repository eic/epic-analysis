// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Ralf Seidl, Christopher Dilks, Sanghwa Park

#ifndef AnalysisEcce_
#define AnalysisEcce_

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "Analysis.h"

class ClustersEE
{
  public:
    ClustersEE() {}
    ClustersEE(double E_, double x_, double y_, double z_, double theta_, double phi_) {}
    virtual ~ClustersEE() {}

    double E;
    double x;
    double y;
    double z;
    double theta;
    double phi;

};

class ParticlesEE
{
  public:
    int pid;
    int charge;
    int mcID;
    TLorentzVector vecPart;
};

class AnalysisEcce : public Analysis
{
  public:
    AnalysisEcce(TString configFileName_="", TString outfilePrefix_="");
    ~AnalysisEcce();

    void Execute() override;

    // select which track source; can be set at the macro level
    UShort_t trackSource; /* all = 0,
                           * inner = 1,
                           * silicon = 2,
                           * ttl = 3
                           */

    ClassDefOverride(AnalysisEcce,1);
};


static const double pimass = 0.13957061;
static const double kmass  = 0.493677;
static const double pmass = 0.938272081;
static const double emass = 0.000511;
static const double mumass = 0.105658376;


#endif
