// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// common constants for analysis

#pragma once

namespace constants {

  // PDG codes
  enum pdg_enum {
    pdgElectron = 11,
    pdgProton   = 2212
  };

  // status codes
  enum status_enum {
    statusFinal = 1,
    statusBeam  = 4
  };

  static const double pimass = 0.13957061;
  static const double kmass  = 0.493677;
  static const double pmass  = 0.938272081;
  static const double emass  = 0.000511;
  static const double mumass = 0.105658376;

};
