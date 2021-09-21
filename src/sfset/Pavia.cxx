#include "Pavia.h"

#include "../interp/Interpolate.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <string>

#define SF_SET_DIR "grids"
#define PAVIA_DIR "pavia"

namespace {

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

std::istream& find_file(std::ifstream& fin, char const* file_name) {
	fin.open(std::string("../share/" SF_SET_DIR "/" PAVIA_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(SF_SET_DIR "/" PAVIA_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(PAVIA_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	throw DataFileNotFoundError(file_name);
}

std::array<Grid<double, 4>, 2> load_grids(char const* file_name) {
	std::ifstream in;
	find_file(in, file_name);
	std::vector<std::array<double, 6> > data;
	std::string line;
	std::getline(in, line);
	while (std::getline(in, line)) {
		std::string blank;
		std::stringstream line_in(line);
		std::array<double, 6> next;
		line_in >> blank >> blank;
		line_in >> next[0] >> next[1] >> next[2] >> next[3];
		line_in >> blank >> blank >> blank >> blank >> blank;
		line_in >> next[4] >> next[5];
		data.push_back(next);
		if (!line_in) {
			throw DataFileParseError(file_name);
		}
	}
	// By default, assume the grids are only accurate to single precision.
	return read_grids<double, 4, 2>(data);
}

}

struct PaviaSfSet::Impl {
	std::array<Grid<double, 4>, 2> grid_pip;
	std::array<Grid<double, 4>, 2> grid_pim;
	CubicView<double, 4> interp_pip_uu;
	CubicView<double, 4> interp_pip_ut_sivers;
	CubicView<double, 4> interp_pim_uu;
	CubicView<double, 4> interp_pim_ut_sivers;
	Impl() :
		grid_pip(load_grids("grid_Pip.txt")),
		grid_pim(load_grids("grid_Pim.txt")),
		interp_pip_uu(grid_pip[0]),
		interp_pip_ut_sivers(grid_pip[1]),
		interp_pim_uu(grid_pim[0]),
		interp_pim_ut_sivers(grid_pim[1]) {
	}
};

PaviaSfSet::PaviaSfSet(PaviaSfSet&& other) noexcept :
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
PaviaSfSet& PaviaSfSet::operator=(PaviaSfSet&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

PaviaSfSet::PaviaSfSet() : SfSet() {
	_impl = new Impl();
}

PaviaSfSet::~PaviaSfSet() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

double PaviaSfSet::F_UUT(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	double Q = std::sqrt(Q_sq);
	double qT = std::sqrt(ph_t_sq) / (z * Q);
	if (h == Hadron::PI_P) {
		return _impl->interp_pip_uu({ Q, x, z, qT });
	} else if (h == Hadron::PI_M) {
		return _impl->interp_pim_uu({ Q, x, z, qT });
	} else {
		return std::numeric_limits<double>::quiet_NaN();
	}
}

double PaviaSfSet::F_UTT_sin_phih_m_phis(Hadron h, double x, double z, double Q_sq, double ph_t_sq) const {
	double Q = std::sqrt(Q_sq);
	double qT = std::sqrt(ph_t_sq) / (z * Q);
	if (h == Hadron::PI_P) {
		return _impl->interp_pip_ut_sivers({ Q, x, z, qT });
	} else if (h == Hadron::PI_M) {
		return _impl->interp_pim_ut_sivers({ Q, x, z, qT });
	} else {
		return std::numeric_limits<double>::quiet_NaN();
	}
}

