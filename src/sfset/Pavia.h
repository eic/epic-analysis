#ifndef Pavia_
#define Pavia_

#include <TMath.h>

#include "../Kinematics.h"
#include "../Weights.h"

class PaviaWeights final : public Weights {
	struct Impl;
	Impl* _impl;

public:
	PaviaWeights();
	PaviaWeights(PaviaWeights const& other) = delete;
	PaviaWeights(PaviaWeights&& other) noexcept;
	PaviaWeights& operator=(PaviaWeights const& other) = delete;
	PaviaWeights& operator=(PaviaWeights&& other) noexcept;
	virtual ~PaviaWeights();

	Double_t GetWeight(const Kinematics& kin) const override;

	// Asymmetries.
	Double_t Sivers(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
	Double_t SiversLower(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
	Double_t SiversUpper(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
};

#endif

