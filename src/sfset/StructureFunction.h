// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#ifndef StructureFunction_
#define StructureFunction_

enum class Hadron {
	PI_0,
	PI_P,
	PI_M,
	K_0,
	K_P,
	K_M,
};

/**
 * Complete set of structure functions bundled together. This abstract class is
 * to be derived for user-provided structure functions. If the structure
 * functions can be factorized into transverse momentum distributions (TMDs),
 * then one of the derived classes TmdSfSet, GaussianTmdSfSet, WwTmdSfSet, or
 * GaussianWwTmdSfSet will be more suitable.
 *
 * This class contains all leading twist and sub-leading twist structure
 * functions.
 */
class SfSet {
public:
	SfSet() { }
	SfSet(SfSet const&) = delete;
	SfSet(SfSet&&) = delete;
	SfSet& operator=(SfSet const&) = delete;
	SfSet& operator=(SfSet&&) = delete;
	virtual ~SfSet() = default;

	virtual double F_UUL(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UU_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UU_cos_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }

	virtual double F_UL_sin_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UL_sin_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }

	virtual double F_UTL_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UT_sin_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UT_sin_3phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UT_sin_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }

	virtual double F_LU_sin_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }

	virtual double F_LL(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_LL_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }

	virtual double F_LT_cos_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_LT_cos_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
	virtual double F_LT_cos_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
        return 0.;
    }
};

#endif

