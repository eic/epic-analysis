#ifndef Allgrid_
#define Allgrid_

#include <TMath.h>

#include "../Kinematics.h"
#include "../Weights.h"

class AllgridWeights final : public Weights {
	struct Impl;
	Impl* _impl;

public:
	AllgridWeights();
	AllgridWeights(AllgridWeights const& other) = delete;
	AllgridWeights(AllgridWeights&& other) noexcept;
	AllgridWeights& operator=(AllgridWeights const& other) = delete;
	AllgridWeights& operator=(AllgridWeights&& other) noexcept;
	virtual ~AllgridWeights();

	Double_t GetWeight(const Kinematics& kin) const override;

	// Asymmetries.
	Double_t ALL(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
	Double_t ALL_lower(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
	Double_t ALL_upper(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const;
};

#endif

