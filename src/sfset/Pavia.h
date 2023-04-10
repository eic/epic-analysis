// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#ifndef Pavia_
#define Pavia_

#include "StructureFunction.h"

class PaviaSfSet final : public SfSet {
	struct Impl;
	Impl* _impl;

public:
	PaviaSfSet();
	PaviaSfSet(PaviaSfSet const& other) = delete;
	PaviaSfSet(PaviaSfSet&& other) noexcept;
	PaviaSfSet& operator=(PaviaSfSet const& other) = delete;
	PaviaSfSet& operator=(PaviaSfSet&& other) noexcept;
	virtual ~PaviaSfSet();

	// Structure functions.
	double F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override;
	double F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const override {
		return 0.;
	}
};

#endif

