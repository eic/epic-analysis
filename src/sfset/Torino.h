#ifndef Torino_
#define Torino_

#include "StructureFunction.h"

class TorinoSfSet final : public SfSet {
	struct Impl;
	Impl* _impl;

public:
	TorinoSfSet();
	TorinoSfSet(TorinoSfSet const& other) = delete;
	TorinoSfSet(TorinoSfSet&& other) noexcept;
	TorinoSfSet& operator=(TorinoSfSet const& other) = delete;
	TorinoSfSet& operator=(TorinoSfSet&& other) noexcept;
	virtual ~TorinoSfSet();

	// Structure functions.
	double F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
};

#endif

