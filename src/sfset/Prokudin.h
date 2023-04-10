// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#ifndef Prokudin_
#define Prokudin_

#include "StructureFunction.h"

/**
 * Use data files from \cite bastami2019ww to calculate structure functions.
 */
class ProkudinSfSet final : public SfSet {
	struct Impl;
	Impl* _impl;

	// Fragmentation functions.
	double D1(Hadron h, unsigned fl, double z, double Q_sq) const;
	double H1perpM1(Hadron h, unsigned fl, double z, double Q_sq) const;

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
	ProkudinSfSet();
	ProkudinSfSet(ProkudinSfSet const& other) = delete;
	ProkudinSfSet(ProkudinSfSet&& other) noexcept;
	ProkudinSfSet& operator=(ProkudinSfSet const& other) = delete;
	ProkudinSfSet& operator=(ProkudinSfSet&& other) noexcept;
	virtual ~ProkudinSfSet();

	// Structure functions.
	double F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UU_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UU_cos_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;

	double F_UL_sin_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UL_sin_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;

	double F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_3phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;

	double F_LL(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_LL_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;

	double F_LT_cos_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_LT_cos_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_LT_cos_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
};

#endif

