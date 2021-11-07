#include "Pavia.h"

#include "../interp/Interpolate.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <string>

#include <TMath.h>

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

std::array<Grid<double, 4>, 3> load_grids(char const* file_name) {
	std::ifstream in;
	find_file(in, file_name);
	std::vector<std::array<double, 7> > data;
	std::string line;
	std::getline(in, line);
	while (std::getline(in, line)) {
		std::string blank;
		std::stringstream line_in(line);
		std::array<double, 7> next;
		line_in >> next[0] >> next[1] >> next[2] >> next[3];
		line_in >> next[4] >> next[5] >> next[6];
		data.push_back(next);
		if (!line_in) {
			throw DataFileParseError(file_name);
		}
	}
	// By default, assume the grids are only accurate to single precision.
	return read_grids<double, 4, 3>(data);
}

}

struct PaviaWeights::Impl {
	std::array<Grid<double, 4>, 3> grid_pip;
	std::array<Grid<double, 4>, 3> grid_pim;
	CubicView<double, 4> interp_pip_sivers;
	CubicView<double, 4> interp_pip_sivers_lower;
	CubicView<double, 4> interp_pip_sivers_upper;
	CubicView<double, 4> interp_pim_sivers;
	CubicView<double, 4> interp_pim_sivers_lower;
	CubicView<double, 4> interp_pim_sivers_upper;
	Impl() :
		grid_pip(load_grids("grid_Pip.txt")),
		grid_pim(load_grids("grid_Pim.txt")),
		interp_pip_sivers(grid_pip[1]),
		interp_pip_sivers_lower(grid_pip[0]),
		interp_pip_sivers_upper(grid_pip[2]),
		interp_pim_sivers(grid_pim[1]),
		interp_pim_sivers_lower(grid_pim[0]),
		interp_pim_sivers_upper(grid_pim[2]) { }
};

PaviaWeights::PaviaWeights(PaviaWeights&& other) noexcept :
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
PaviaWeights& PaviaWeights::operator=(PaviaWeights&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

PaviaWeights::PaviaWeights() : Weights() {
	_impl = new Impl();
}

PaviaWeights::~PaviaWeights() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

Double_t PaviaWeights::GetWeight(const Kinematics& kin) const {
	Double_t sivers = Sivers(kin.hadPID, kin.x, kin.z, kin.Q2, kin.pT);
	if (!TMath::Finite(sivers) || TMath::Abs(sivers) > 1.) {
		sivers = 0.;
	}
	return 1. + kin.polT * kin.tSpin
		* sivers * TMath::Sin(kin.phiH - kin.phiS);
}

Double_t PaviaWeights::Sivers(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const {
	Double_t Q = TMath::Sqrt(Q2);
	Double_t qTtoQ = pT / (z * Q);
	if (hadPID == 211) {
		return _impl->interp_pip_sivers({ Q, x, z, qTtoQ });
	} else if (hadPID == -211) {
		return _impl->interp_pim_sivers({ Q, x, z, qTtoQ });
	} else {
		return 0.;
	}
}

Double_t PaviaWeights::SiversLower(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const {
	Double_t Q = TMath::Sqrt(Q2);
	Double_t qTtoQ = pT / (z * Q);
	if (hadPID == 211) {
		return _impl->interp_pip_sivers_lower({ Q, x, z, qTtoQ });
	} else if (hadPID == -211) {
		return _impl->interp_pim_sivers_lower({ Q, x, z, qTtoQ });
	} else {
		return 0.;
	}
}

Double_t PaviaWeights::SiversUpper(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const {
	Double_t Q = TMath::Sqrt(Q2);
	Double_t qTtoQ = pT / (z * Q);
	if (hadPID == 211) {
		return _impl->interp_pip_sivers_upper({ Q, x, z, qTtoQ });
	} else if (hadPID == -211) {
		return _impl->interp_pim_sivers_upper({ Q, x, z, qTtoQ });
	} else {
		return 0.;
	}
}

