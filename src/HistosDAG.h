// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Christopher Dilks

// Implement Histos payload for Adage
#pragma once

#include <adage/Adage.h>
#include "Histos.h"

class HistosDAG : public Adage<Histos>
{
  public:
    HistosDAG() : Adage<Histos>("histos") {};
    ~HistosDAG() {};

    // build the DAG from specified bin scheme
    void Build(std::map<TString,BinSet*> binSchemes) override;

  ClassDefOverride(HistosDAG,1);
};
