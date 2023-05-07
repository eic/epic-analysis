// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Sanghwa Park, Christopher Dilks

#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "Analysis.h"

class AnalysisAthena : public Analysis
{
  public:
    AnalysisAthena(TString configFileName_="", TString outfilePrefix_="");
    ~AnalysisAthena();

    void Execute() override;

    ClassDefOverride(AnalysisAthena,1);
};
