// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Christopher Dilks

#pragma once

#include <vector>
#include "Analysis.h"
#include "Kinematics.h"

class Analysis;

class Dihadron {
  public:
    Dihadron() {};
    Dihadron(TLorentzVector vecHad1, TLorentzVector vecHad2, Kinematics *K);
    ~Dihadron() {};
    Double_t Mh, MX, Z, Ph, PhPerp, Theta, PhiH, PhiR, PhiS, Zeta;
    Double_t hadZ[2];
  private:
  ClassDef(Dihadron,1);
};

class DihadronSet {
  public:
    DihadronSet() : debug(false), dihadronFinalStateID("") {};
    ~DihadronSet() {};
    void IncludeHadron(TString hadName);
    void AddHadron(Analysis *A);
    void SetFinalStateID(TString state) { dihadronFinalStateID = state; }
    void CalculateKinematics(Analysis *A, Double_t wgt=1.0);
  private:
    std::vector<TString> includedHadrons;
    std::map<TString,std::vector<TLorentzVector>> hadSetRec, hadSetGen;
    std::vector<Dihadron> dihListRec, dihListGen;
    TString dihadronFinalStateID;
    bool debug;
  ClassDef(DihadronSet,1);
};

