// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Sanghwa Park, Christopher Dilks

#pragma once

#include <TLorentzVector.h>

class Particles {
  public:
    int pid    = 0;
    int charge = 0;
    int mcID   = -1;
    TLorentzVector vecPart;
};

/*
class Clusters {
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
*/

