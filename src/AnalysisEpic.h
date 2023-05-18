// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Gregory Matousek, Christopher Dilks

#pragma once

// ROOT
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// epic-analysis
#include "Analysis.h"

class AnalysisEpic : public Analysis
{
  public:
    AnalysisEpic(TString infileName_="", TString outfilePrefix_="");
    ~AnalysisEpic();

    void Execute() override;

    ClassDefOverride(AnalysisEpic,1);
};
