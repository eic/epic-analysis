// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Ralf Seidl, Christopher Dilks, Sanghwa Park

#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "Analysis.h"

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
