#include "Allgrid.h"

#include "../interp/Interpolate.h"

#include <fstream>
#include <stdexcept>
#include <sstream>
#include <string>

#include <TMath.h>

#define SF_SET_DIR "grids"
#define ALLGRID_DIR "allgrid"

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
	fin.open(std::string("../share/" SF_SET_DIR "/" ALLGRID_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(SF_SET_DIR "/" ALLGRID_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(std::string(ALLGRID_DIR "/") + file_name);
	if (fin) {
		return fin;
	}
	fin.open(file_name);
	if (fin) {
		return fin;
	}
	throw DataFileNotFoundError(file_name);
}

std::array<Grid<double, 3>, 1> load_grids(char const* file_name) {
	std::ifstream in;
	find_file(in, file_name);
	std::vector<std::array<double, 4> > data;
	std::string line;
	while (std::getline(in, line)) {
		std::string blank;
		std::stringstream line_in(line);
		std::array<double, 4> next;
		line_in >> next[0] >> next[1] >> next[2];
		line_in >> next[3] >> blank;
		data.push_back(next);
		if (!line_in) {
			throw DataFileParseError(file_name);
		}
	}
	// By default, assume the grids are only accurate to single precision.
	return read_grids<double, 3, 1>(data);
}

}

struct AllgridWeights::Impl {
	std::array<Grid<double, 3>, 1> grid_pip;
	CubicView<double, 3> interp_pip;
	Impl() :
		grid_pip(load_grids("grid_Pip.txt")),
		interp_pip(grid_pip[0]) { }
};

AllgridWeights::AllgridWeights(AllgridWeights&& other) noexcept :
		_impl(nullptr) {
	std::swap(_impl, other._impl);
}
AllgridWeights& AllgridWeights::operator=(AllgridWeights&& other) noexcept {
	std::swap(_impl, other._impl);
	return *this;
}

AllgridWeights::AllgridWeights() : Weights() {
	_impl = new Impl();
}

AllgridWeights::~AllgridWeights() {
	if (_impl != nullptr) {
		delete _impl;
	}
}

Double_t AllgridWeights::GetWeight(const Kinematics& kin) const {
	Double_t all = ALL(kin.hadPID, kin.x, kin.z, kin.Q2, kin.pT);
	if (!TMath::Finite(all) || TMath::Abs(all) > 1.) {
		all = 0.;
	}
	return kin.polL * kin.polBeam * kin.lSpin * kin.depolP2 * all;
}

Double_t AllgridWeights::ALL(Int_t hadPID, Double_t x, Double_t z, Double_t Q2, Double_t pT) const {
	if (hadPID == 211) {
		return _impl->interp_pip({ Q2, x, z });
	} else {
		return 0.;
	}
}

