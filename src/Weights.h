#ifndef Weights_
#define Weights_

#include <memory>
#include <initializer_list>
#include <vector>

#include "Kinematics.h"

// Interface for all weight providers.
class Weights : public TNamed
{
  public:
    Weights() { }
    // Return how an event should be weighted based on its kinematics.
    virtual Double_t GetWeight(const Kinematics&) const = 0;
  private:
  ClassDef(Weights,1);
};

// Provide a constant weight to all events.
class WeightsUniform : public Weights
{
  public:
    WeightsUniform(Double_t weight = 1.) : weight(weight) { }
    Double_t GetWeight(const Kinematics&) const override {
        return weight;
    }
  private:
    Double_t weight;
  ClassDefOverride(WeightsUniform,1);
};

// Inject a Sivers asymmetry.
class WeightsSivers : public Weights
{
  public:
    Double_t GetWeight(const Kinematics& kin) const override;
    virtual Double_t Asymmetry(Int_t h, Double_t x, Double_t z, Double_t Q2, Double_t pt) const = 0;
  private:
  ClassDefOverride(WeightsSivers,1);
};

// Inject a Collins asymmetry.
class WeightsCollins : public Weights
{
  public:
    Double_t GetWeight(const Kinematics& kin) const override;
    virtual Double_t Asymmetry(Int_t h, Double_t x, Double_t z, Double_t Q2, Double_t pt) const = 0;

  private:
  ClassDefOverride(WeightsCollins,1);
};

// Inject a ALL asymmetry.
class WeightsALL : public Weights
{
  public:
    Double_t GetWeight(const Kinematics& kin) const override;
    virtual Double_t Asymmetry(Int_t h, Double_t x, Double_t z, Double_t Q2, Double_t pt) const = 0;
  private:
  ClassDefOverride(WeightsALL,1);
};

// Product of weights. Useful for products of acceptances.
class WeightsProduct : public Weights
{
  public:
    WeightsProduct(std::initializer_list<Weights const*> weights);

    Double_t GetWeight(const Kinematics& kin) const override;
    WeightsProduct& Multiply(Weights const* rhs);

  private:
    std::vector<Weights const*> weights;
  ClassDefOverride(WeightsProduct,1);
};

// Sum of weights. Useful for adding asymmetries together.
class WeightsSum : public Weights
{
  public:
    WeightsSum(std::initializer_list<Weights const*> weights);

    Double_t GetWeight(const Kinematics& kin) const override;
    WeightsSum& Add(Weights const* rhs);

  private:
    std::vector<Weights const*> weights;
  ClassDefOverride(WeightsSum,1);
};

#endif
