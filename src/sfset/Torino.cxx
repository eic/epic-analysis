#include "Torino.h"

#include <cmath>
#include <string>

#include <stfunctions.h>
#include <parameters.h>

struct TorinoSfSet::Impl {
	TMD torino;
	Impl() { }
};

TorinoSfSet::TorinoSfSet(TorinoSfSet&& other) noexcept :
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
TorinoSfSet& TorinoSfSet::operator=(TorinoSfSet&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

TorinoSfSet::TorinoSfSet() : SfSet() {
	_impl = new Impl();
}

TorinoSfSet::~TorinoSfSet() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

char const* hadron_name(Hadron h) {
	switch (h) {
	case Hadron::PI_P:
		return "pi+";
	case Hadron::PI_M:
		return "pi-";
	}
	return "<error>";
}

double TorinoSfSet::F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	// The Torino calculation takes S for some reason, but doesn't depend on the
	// what value of S is provided. So just make up some semi-reasonable value
	// for it.
	double S = 0.5 * Q_sq / x;
	return _impl->torino.FUU("proton", hadron_name(h), S, x, z, Q_sq, std::sqrt(ph_t_sq));
}

double TorinoSfSet::F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	double S = 0.5 * Q_sq / x;
	return _impl->torino.FUTSivers("proton", hadron_name(h), S, x, z, Q_sq, std::sqrt(ph_t_sq));
}

double TorinoSfSet::F_UT_sin_phih_p_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	double S = 0.5 * Q_sq / x;
	return _impl->torino.FUTCollins("proton", hadron_name(h), S, x, z, Q_sq, std::sqrt(ph_t_sq));
}

