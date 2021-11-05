#ifndef Prokudin_
#define Prokudin_

#include "../Weights.h"

/**
 * Use data files from \cite bastami2019ww to calculate structure functions.
 */
class ProkudinWeights final : public Weights {
	struct Impl;
	Impl* _impl;

	// Fragmentation functions.
	double D1(Int_t h, unsigned fl, double z, double Q_sq) const;
	double H1perpM1(Int_t h, unsigned fl, double z, double Q_sq) const;

	// Transverse momentum distributions.
	double xf1(unsigned fl, double x, double Q_sq) const;
	double xf1TperpM1(unsigned fl, double x, double Q_sq) const;
	double xg1(unsigned fl, double x, double Q_sq) const;
	double xgT(unsigned fl, double x, double Q_sq) const;
	double xh1(unsigned fl, double x, double Q_sq) const;
	double xh1M1(unsigned fl, double x, double Q_sq) const;
	double xh1LperpM1(unsigned fl, double x, double Q_sq) const;
	double xh1TperpM2(unsigned fl, double x, double Q_sq) const;
	double xh1perpM1(unsigned fl, double x, double Q_sq) const;

public:
	ProkudinWeights();
	ProkudinWeights(ProkudinWeights const& other) = delete;
	ProkudinWeights(ProkudinWeights&& other) noexcept;
	ProkudinWeights& operator=(ProkudinWeights const& other) = delete;
	ProkudinWeights& operator=(ProkudinWeights&& other) noexcept;
	virtual ~ProkudinWeights();

	double GetWeight(const Kinematics& kin) const override;

	// Structure functions.
	double F_UUT(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UU_cos_phih(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UU_cos_2phih(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;

	double F_UL_sin_phih(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UL_sin_2phih(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;

	double F_UTT_sin_phih_m_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UT_sin_2phih_m_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UT_sin_3phih_m_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UT_sin_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_UT_sin_phih_p_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;

	double F_LL(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_LL_cos_phih(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;

	double F_LT_cos_phih_m_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_LT_cos_2phih_m_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
	double F_LT_cos_phis(Int_t h, double x, double z, double Q_sq, double ph_t_sq) const;
};

#endif

