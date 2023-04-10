// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#include "Prokudin.h"

#include "../interp/Interpolate.h"

#include <mstwpdf.h>

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#define SF_SET_DIR "grids"
#define PROKUDIN_DIR "prokudin"

// The following parameters come from appendix A in [2].

namespace {

unsigned const NUM_FLAVORS = 6;

// We only have data files for proton structure function.
double const M  = 0.93827208816;
double const mh = 0.13957039;
double const PI = 3.1415926535897932384626433;
double const E  = 2.7182818284590452353602875;

// Parameters for PDF.
double const F1_MEAN_K_PERP_SQ = 0.25;
double const G1_MEAN_K_PERP_SQ = 0.76 * F1_MEAN_K_PERP_SQ;
double const D1_MEAN_P_PERP_SQ = 0.2;

// Parameters for transversity.
double const H1_ALPHA = 1.11;
double const H1_BETA = 3.64;
double const H1_N[6] = {
	0.46, -1.00, 0.,
	0., 0., 0.,
};
double const H1_MEAN_K_PERP_SQ = 0.25;
// Reference [2] is not consistent about the choice of `mean_hT`, so here we
// choose something in between the two options.
double const HT_MEAN_K_PERP_SQ = 0.15;

// Parameters for Collins.
double const COLLINS_M_SQ = 1.50;
double const COLLINS_M = std::sqrt(COLLINS_M_SQ);
double const COLLINS_GAMMA = 1.06;
double const COLLINS_DELTA = 0.07;
double const COLLINS_N_FAV = 0.49;
double const COLLINS_N_DISFAV = -1.;
double const COLLINS_MEAN_P_PERP_SQ = (D1_MEAN_P_PERP_SQ * COLLINS_M_SQ)
	/ (D1_MEAN_P_PERP_SQ + COLLINS_M_SQ);

// Parameters for Sivers.
double const SIVERS_M_1_SQ = 0.19;
double const SIVERS_M_1 = std::sqrt(SIVERS_M_1_SQ);
double const SIVERS_ALPHA[6] = {
	0.35, 0.44, 0.,
	0., 0., 0.,
};
double const SIVERS_BETA[6] = {
	2.6, 0.90, 0.,
	0., 0., 0.,
};
double const SIVERS_N[6] = {
	0.40, -0.97, 0.,
	0., 0., 0.,
};
double const SIVERS_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * SIVERS_M_1_SQ)
	/ (F1_MEAN_K_PERP_SQ + SIVERS_M_1_SQ);

// Parameters for Boer-Mulders.
double const BM_M_1_SQ = 0.34;
double const BM_M_1 = std::sqrt(BM_M_1_SQ);
double const BM_ALPHA[6] = {
	0.73, 1.08, 0.79,
	0.79, 0.79, 0.79,
};
double const BM_BETA = 3.46;
double const BM_A[6] = {
	0.35, -0.90, 0.24,
	0.04, -0.40, -1.,
};
double const BM_LAMBDA[6] = {
	2.1, -1.111, 0.,
	0., 0., 0.,
};
double const BM_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * BM_M_1_SQ)
	/ (F1_MEAN_K_PERP_SQ + BM_M_1_SQ);

// Parameters for pretzelosity.
double const PRETZ_M_TT_SQ = 0.18;
double const PRETZ_M_TT = std::sqrt(PRETZ_M_TT_SQ);
double const PRETZ_ALPHA = 2.5;
double const PRETZ_BETA = 2.;
double const PRETZ_N[6] = {
	1., -1., 0.,
	0., 0., 0.,
};

double sq(double x) {
	return x * x;
}
double const PRETZ_MEAN_K_PERP_SQ = (F1_MEAN_K_PERP_SQ * PRETZ_M_TT_SQ)
	/ (F1_MEAN_K_PERP_SQ + PRETZ_M_TT_SQ);

double G(double ph_t_sq, double l) {
	// Equation [2.5.2].
	return std::exp(-ph_t_sq / l) / (PI * l);
}
double lambda(double z, double mean_kperp_sq, double mean_pperp_sq) {
	return sq(z) * mean_kperp_sq + mean_pperp_sq;
}

struct DataFileNotFoundError : public std::runtime_error {
	DataFileNotFoundError(char const* name) :
		std::runtime_error(std::string("Couldn't load data file '") + name + "'") {
	}
};
struct DataFileParseError : public std::runtime_error {
	DataFileParseError(char const* name) :
		std::runtime_error(std::string("Couldn't parse data file '") + name + "'") {
	}
};

// Finds a grid file.
std::istream& find_file(std::ifstream& fin, char const* file_name) {
	fin.open(std::string("../share/" SF_SET_DIR "/" PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(SF_SET_DIR "/" PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(PROKUDIN_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	throw DataFileNotFoundError(file_name);
}

// Convenience method for loading grid data from a file. Automatically searches
// in several possible directories for the file.
template<typename T, std::size_t N, std::size_t K>
std::array<Grid<T, N>, K> load_grids(char const* file_name) {
	std::ifstream in;
	find_file(in, file_name);
	std::vector<std::array<T, N + K> > data;
	std::string line;
	while (std::getline(in, line)) {
		std::stringstream line_in(line);
		std::array<T, N + K> next;
		for (std::size_t idx = 0; idx < N + K; ++idx) {
			line_in >> next[idx];
		}
		data.push_back(next);
		if (!line_in) {
			throw DataFileParseError(file_name);
		}
	}
	// By default, assume the grids are only accurate to single precision.
	return read_grids<T, N, K>(data);
}

double charge(unsigned fl) {
	switch (fl) {
	// Up.
	case 0:
		return 2./3.;
	// Down.
	case 1:
	// Strange.
	case 2:
		return -1./3.;
	// Up bar.
	case 3:
		return -2./3.;
	// Down bar.
	case 4:
	// Strange bar.
	case 5:
		return 1./3.;
	default:
		return 0.;
	}
}

struct ProkudinImpl {
	// PDF are calculated with the MSTWPDF library.
	std::ifstream file_pdf;
	mstw::c_mstwpdf pdf;

	// Load data files from WWSIDIS repository for the TMDs and FFs.
	std::array<Grid<double, 2>, 6> data_D1_pi_plus;
	std::array<Grid<double, 2>, 6> data_D1_pi_minus;
	std::array<Grid<double, 2>, 6> data_g1;
	std::array<Grid<double, 2>, 6> data_xgT;
	std::array<Grid<double, 2>, 2> data_xh1LperpM1;
	// Soffer bound.
	std::array<Grid<double, 2>, 6> data_sb;

	// Fragmentation functions.
	CubicView<double, 2> interp_D1_pi_plus[6];
	CubicView<double, 2> interp_D1_pi_minus[6];

	// Transverse momentum distributions.
	CubicView<double, 2> interp_g1[6];
	CubicView<double, 2> interp_xgT[6];
	CubicView<double, 2> interp_xh1LperpM1[2];
	CubicView<double, 2> interp_sb[6];

	ProkudinImpl() :
			file_pdf(),
			pdf(find_file(file_pdf, "mstw2008lo.00.dat"), false, true),
			data_D1_pi_plus(
				load_grids<double, 2, 6>("fragmentationpiplus.dat")),
			data_D1_pi_minus(
				load_grids<double, 2, 6>("fragmentationpiminus.dat")),
			data_g1(load_grids<double, 2, 6>("g1.dat")),
			data_xgT(load_grids<double, 2, 6>("gT_u_d_ubar_dbar_s_sbar.dat")),
			data_xh1LperpM1(load_grids<double, 2, 2>("xh1Lperp_u_d.dat")),
			data_sb(load_grids<double, 2, 6>("SofferBound.dat")),
			interp_D1_pi_plus{
				CubicView<double, 2>(data_D1_pi_plus[0]),
				CubicView<double, 2>(data_D1_pi_plus[1]),
				CubicView<double, 2>(data_D1_pi_plus[2]),
				CubicView<double, 2>(data_D1_pi_plus[3]),
				CubicView<double, 2>(data_D1_pi_plus[4]),
				CubicView<double, 2>(data_D1_pi_plus[5]),
			},
			interp_D1_pi_minus{
				CubicView<double, 2>(data_D1_pi_minus[0]),
				CubicView<double, 2>(data_D1_pi_minus[1]),
				CubicView<double, 2>(data_D1_pi_minus[2]),
				CubicView<double, 2>(data_D1_pi_minus[3]),
				CubicView<double, 2>(data_D1_pi_minus[4]),
				CubicView<double, 2>(data_D1_pi_minus[5]),
			},
			interp_g1{
				CubicView<double, 2>(data_g1[0]),
				CubicView<double, 2>(data_g1[1]),
				CubicView<double, 2>(data_g1[2]),
				CubicView<double, 2>(data_g1[3]),
				CubicView<double, 2>(data_g1[4]),
				CubicView<double, 2>(data_g1[5]),
			},
			interp_xgT{
				CubicView<double, 2>(data_xgT[0]),
				CubicView<double, 2>(data_xgT[1]),
				CubicView<double, 2>(data_xgT[4]),
				CubicView<double, 2>(data_xgT[2]),
				CubicView<double, 2>(data_xgT[3]),
				CubicView<double, 2>(data_xgT[5]),
			},
			interp_xh1LperpM1{
				CubicView<double, 2>(data_xh1LperpM1[0]),
				CubicView<double, 2>(data_xh1LperpM1[1]),
			},
			interp_sb{
				CubicView<double, 2>(data_sb[0]),
				CubicView<double, 2>(data_sb[1]),
				CubicView<double, 2>(data_sb[2]),
				CubicView<double, 2>(data_sb[3]),
				CubicView<double, 2>(data_sb[4]),
				CubicView<double, 2>(data_sb[5]),
			} {
		file_pdf.close();
	}
};

}

struct ProkudinSfSet::Impl {
	ProkudinImpl impl;
};

ProkudinSfSet::ProkudinSfSet(ProkudinSfSet&& other) noexcept :
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
ProkudinSfSet& ProkudinSfSet::operator=(ProkudinSfSet&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

ProkudinSfSet::ProkudinSfSet() : SfSet() {
	_impl = new Impl();
}

ProkudinSfSet::~ProkudinSfSet() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

double ProkudinSfSet::F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.1a].
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xf1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	double l = lambda(z, F1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UU_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.9a].
	double Q = std::sqrt(Q_sq);
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	// Uses a WW-type approximation to rewrite in terms of `xf1`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xf1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	double l = lambda(z, F1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*F1_MEAN_K_PERP_SQ/Q*ph_t*(z/l)*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UU_cos_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.9a].
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1perpM1(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	double l = lambda(z, BM_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 4.*M*mh*ph_t_sq*sq(z/l)*G(ph_t_sq, l)*result;
}

double ProkudinSfSet::F_UL_sin_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.6a].
	double Q = std::sqrt(Q_sq);
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	// Use WW-type approximation to rewrite in terms of `xh1LperpM1`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1LperpM1(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	// Approximate width with `H1_MEAN_K_PERP_SQ`.
	double l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return -8.*M*mh*z*ph_t/(Q*l)*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UL_sin_2phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.6.2a].
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1LperpM1(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	// Approximate width with `H1_MEAN_K_PERP_SQ`.
	double l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 4.*M*mh*ph_t_sq*sq(z/l)*G(ph_t_sq, l)*result;
}

double ProkudinSfSet::F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.7a].
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xf1TperpM1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	double l = lambda(z, SIVERS_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*M*z*ph_t/l*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UT_sin_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.8a].
	double Q = std::sqrt(Q_sq);
	double result_1 = 0.;
	// Use WW-type approximation to rewrite in terms of `xf1TperpM1`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result_1 += sq(charge(fl))*xf1TperpM1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	double result_2 = 0.;
	// Use WW-type approximation to rewrite in terms of `h1TperpM2`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result_2 += sq(charge(fl))*xh1TperpM2(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	// Approximate width with `SIVERS_MEAN_K_PERP_SQ`.
	double l_1 = lambda(z, SIVERS_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	// Approximate width with `PRETZ_MEAN_K_PERP_SQ`.
	double l_2 = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	// The paragraph following [2.7.8a] has a mistake in the WW-type
	// approximation linking `h1TM1 + h1TperpM1` with `h1TperpM2`, due to a
	// missing factor of 2.
	return 2.*M*ph_t_sq/Q*(
		SIVERS_MEAN_K_PERP_SQ*sq(z/l_1)*G(ph_t_sq, l_1)*result_1
		- 2.*M*mh*sq(z/l_2)*G(ph_t_sq, l_2)*result_2);
}
double ProkudinSfSet::F_UT_sin_3phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.10a].
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1TperpM2(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	double l = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 2.*sq(M)*mh*std::pow(z*ph_t/l, 3)*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UT_sin_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.7a].
	double Q = std::sqrt(Q_sq);
	double result = 0.;
	// WW-type approximation used here (see [2] for details).
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1M1(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	double l = lambda(z, PRETZ_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 8.*sq(M)*mh*sq(z)/(Q*l)*(1. - ph_t_sq/l)*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.8a].
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xh1(fl, x, Q_sq)*H1perpM1(h, fl, z, Q_sq);
	}
	double l = lambda(z, H1_MEAN_K_PERP_SQ, COLLINS_MEAN_P_PERP_SQ);
	return 2.*mh*z*ph_t/l*G(ph_t_sq, l)*result;
}

double ProkudinSfSet::F_LL(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.5.5a].
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xg1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	double l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_LL_cos_phih(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.5a].
	double Q = std::sqrt(Q_sq);
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	// Uses a WW-type approximation to rewrite in terms of `xg1`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xg1(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	double l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*G1_MEAN_K_PERP_SQ*z*ph_t/(Q*l)*G(ph_t_sq, l)*result;
}

double ProkudinSfSet::F_LT_cos_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.6.1a].
	double ph_t = std::sqrt(ph_t_sq);
	double result = 0.;
	// Uses a WW-type approximation to rewrite in terms of `xgT`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xgT(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	double l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return 2.*M*x*z*ph_t/l*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_LT_cos_2phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.4a].
	double Q = std::sqrt(Q_sq);
	double result = 0.;
	// Uses a WW-type approximation to rewrite in terms of `xgT`.
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xgT(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	double l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*G1_MEAN_K_PERP_SQ*M*x*ph_t_sq*sq(z/l)/Q*G(ph_t_sq, l)*result;
}
double ProkudinSfSet::F_LT_cos_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// Equation [2.7.2a].
	double Q = std::sqrt(Q_sq);
	double result = 0.;
	for (unsigned fl = 0; fl < NUM_FLAVORS; ++fl) {
		result += sq(charge(fl))*xgT(fl, x, Q_sq)*D1(h, fl, z, Q_sq);
	}
	// Approximate width with `G1_MEAN_K_PERP_SQ`.
	double l = lambda(z, G1_MEAN_K_PERP_SQ, D1_MEAN_P_PERP_SQ);
	return -2.*M*x/Q*G(ph_t_sq, l)*result;
}

// Fragmentation functions.
double ProkudinSfSet::D1(Hadron h, unsigned fl, double z, double Q_sq) const {
	switch (h) {
	case Hadron::PI_P:
		return _impl->impl.interp_D1_pi_plus[fl]({ z, Q_sq });
	case Hadron::PI_M:
		return _impl->impl.interp_D1_pi_minus[fl]({ z, Q_sq });
	default:
		return 0.;
	}
}
double ProkudinSfSet::H1perpM1(Hadron h, unsigned fl, double z, double Q_sq) const {
	double collins_coeff = 0.;
	if (h == Hadron::PI_P) {
		if (fl == 0 || fl == 4) {
			collins_coeff = COLLINS_N_FAV;
		} else if (fl == 1 || fl == 3) {
			collins_coeff = COLLINS_N_DISFAV;
		}
	} else if (h == Hadron::PI_M) {
		if (fl == 1 || fl == 3) {
			collins_coeff = COLLINS_N_FAV;
		} else if (fl == 0 || fl == 4) {
			collins_coeff = COLLINS_N_DISFAV;
		}
	} else {
		return 0.;
	}
	return std::sqrt(E/2.)/(z*mh*COLLINS_M)
		*sq(COLLINS_MEAN_P_PERP_SQ)/D1_MEAN_P_PERP_SQ
		*collins_coeff
		*std::pow(z, COLLINS_GAMMA)*std::pow(1. - z, COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA + COLLINS_DELTA, COLLINS_GAMMA + COLLINS_DELTA)
		*std::pow(COLLINS_GAMMA, -COLLINS_GAMMA)
		*std::pow(COLLINS_DELTA, -COLLINS_DELTA)
		*D1(h, fl, z, Q_sq);
}

// Parton distribution functions.
double ProkudinSfSet::xf1(unsigned fl, double x, double Q_sq) const {
	double Q = std::sqrt(Q_sq);
	switch (fl) {
	case 0:
		return _impl->impl.pdf.parton(8, x, Q) + _impl->impl.pdf.parton(-2, x, Q);
	case 1:
		return _impl->impl.pdf.parton(7, x, Q) + _impl->impl.pdf.parton(-1, x, Q);
	case 2:
		return _impl->impl.pdf.parton(3, x, Q);
	case 3:
		return _impl->impl.pdf.parton(-2, x, Q);
	case 4:
		return _impl->impl.pdf.parton(-1, x, Q);
	case 5:
		return _impl->impl.pdf.parton(-3, x, Q);
	default:
		return 0.;
	}
}

// Transverse momentum distributions.
double ProkudinSfSet::xf1TperpM1(unsigned fl, double x, double Q_sq) const {
	// Equation [2.A.4].
	return -std::sqrt(E/2.)/(M*SIVERS_M_1)
		*sq(SIVERS_MEAN_K_PERP_SQ)/F1_MEAN_K_PERP_SQ
		*SIVERS_N[fl]
		*std::pow(x, SIVERS_ALPHA[fl])*std::pow(1. - x, SIVERS_BETA[fl])
		*std::pow(SIVERS_ALPHA[fl] + SIVERS_BETA[fl],
			SIVERS_ALPHA[fl] + SIVERS_BETA[fl])
		*std::pow(SIVERS_ALPHA[fl], -SIVERS_ALPHA[fl])
		*std::pow(SIVERS_BETA[fl], -SIVERS_BETA[fl])
		*xf1(fl, x, Q_sq);
}
double ProkudinSfSet::xg1(unsigned fl, double x, double Q_sq) const {
	return x*_impl->impl.interp_g1[fl]({ x, Q_sq });
}
double ProkudinSfSet::xgT(unsigned fl, double x, double Q_sq) const {
	return _impl->impl.interp_xgT[fl]({ x, Q_sq });
}
double ProkudinSfSet::xh1(unsigned fl, double x, double Q_sq) const {
	// Use the Soffer bound to get an upper limit on transversity (Equation
	// [2.A.7]).
	return x*H1_N[fl]
		*std::pow(x, H1_ALPHA)*std::pow(1. - x, H1_BETA)
		*std::pow(H1_ALPHA + H1_BETA, H1_ALPHA + H1_BETA)
		*std::pow(H1_ALPHA, -H1_ALPHA)
		*std::pow(H1_BETA, -H1_BETA)
		*_impl->impl.interp_sb[fl]({ x, Q_sq });
}
double ProkudinSfSet::xh1M1(unsigned fl, double x, double Q_sq) const {
	return H1_MEAN_K_PERP_SQ/(2.*sq(M))*xh1(fl, x, Q_sq);
}
double ProkudinSfSet::xh1LperpM1(unsigned fl, double x, double Q_sq) const {
	// Data only exists for up and down quarks.
	if (!(fl == 0 || fl == 1)) {
		return 0.;
	} else {
		return _impl->impl.interp_xh1LperpM1[fl]({ x, Q_sq });
	}
}
double ProkudinSfSet::xh1TperpM2(unsigned fl, double x, double Q_sq) const {
	// Equation [2.A.24].
	return E/(2.*sq(M)*PRETZ_M_TT_SQ)
		*std::pow(PRETZ_MEAN_K_PERP_SQ, 3)/F1_MEAN_K_PERP_SQ
		*PRETZ_N[fl]
		*std::pow(x, PRETZ_ALPHA)*std::pow(1. - x, PRETZ_BETA)
		*std::pow(PRETZ_ALPHA + PRETZ_BETA, PRETZ_ALPHA + PRETZ_BETA)
		*std::pow(PRETZ_ALPHA, -PRETZ_ALPHA)
		*std::pow(PRETZ_BETA, -PRETZ_BETA)
		*(xf1(fl, x, Q_sq) - xg1(fl, x, Q_sq));
}
double ProkudinSfSet::xh1perpM1(unsigned fl, double x, double Q_sq) const {
	// Equation [2.A.18].
	return -std::sqrt(E/2.)/(M*BM_M_1)
		*sq(BM_MEAN_K_PERP_SQ)/F1_MEAN_K_PERP_SQ
		*BM_LAMBDA[fl] * BM_A[fl]
		*std::pow(x, BM_ALPHA[fl])*std::pow(1. - x, BM_BETA)
		*std::pow(BM_ALPHA[fl] + BM_BETA, BM_ALPHA[fl] + BM_BETA)
		*std::pow(BM_ALPHA[fl], -BM_ALPHA[fl])
		*std::pow(BM_BETA, -BM_BETA)
		*xf1(fl, x, Q_sq);
}

