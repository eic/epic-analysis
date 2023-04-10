// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#include "Weights.h"

using std::vector;

ClassImp(Weights)
ClassImp(WeightsUniform)
ClassImp(WeightsSivers)
ClassImp(WeightsCollins)
ClassImp(WeightsProduct)
ClassImp(WeightsSum)

Double_t WeightsSivers::GetWeight(const Kinematics& kin) const {
    return kin.polT * kin.tSpin
        * this->Asymmetry(kin.x, kin.z, kin.Q2, kin.pT)
        * TMath::Sin(kin.phiH - kin.phiS);
}

Double_t WeightsCollins::GetWeight(const Kinematics& kin) const {
    return kin.polT * kin.tSpin * kin.depolP1
        * this->Asymmetry(kin.x, kin.z, kin.Q2, kin.pT)
        * TMath::Sin(kin.phiH + kin.phiS);
}

WeightsProduct::WeightsProduct(std::initializer_list<Weights const*> init) {
    for (Weights const* weight : init) {
        weights.push_back(weight);
    }
}

Double_t WeightsProduct::GetWeight(const Kinematics& kin) const {
    Double_t product = 1.;
    for (auto weight_ptr : weights) {
        product *= weight_ptr->GetWeight(kin);
    }
    return product;
}

WeightsProduct& WeightsProduct::Multiply(Weights const* rhs) {
    WeightsProduct const* rhs_product = dynamic_cast<WeightsProduct const*>(rhs);
    if (rhs_product != nullptr) {
        weights.insert(
            weights.end(),
            rhs_product->weights.begin(),
            rhs_product->weights.end());
    } else {
        weights.push_back(rhs);
    }
    return *this;
}

WeightsSum::WeightsSum(std::initializer_list<Weights const*> init) {
    for (Weights const* weight : init) {
        weights.push_back(weight);
    }
}

Double_t WeightsSum::GetWeight(const Kinematics& kin) const {
    Double_t sum = 0.;
    for (auto weight_ptr : weights) {
        sum += weight_ptr->GetWeight(kin);
    }
    return sum;
}

WeightsSum& WeightsSum::Add(Weights const* rhs) {
    WeightsSum const* rhs_sum = dynamic_cast<WeightsSum const*>(rhs);
    if (rhs_sum != nullptr) {
        weights.insert(
            weights.end(),
            rhs_sum->weights.begin(),
            rhs_sum->weights.end());
    } else {
        weights.push_back(rhs);
    }
    return *this;
}

