// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#include "Hist4D.h"

ClassImp(Hist4D);

struct DifferentNumberOfBins : public std::runtime_error {
	DifferentNumberOfBins() :
		std::runtime_error("histograms have different numbers of bins") { }
};
struct DifferentAxisLimits : public std::runtime_error {
	DifferentAxisLimits() :
		std::runtime_error("histograms have different axis limits") { }
};
struct DifferentBinLimits : public std::runtime_error {
	DifferentBinLimits() :
		std::runtime_error("histograms have different bin limits") { }
};
struct DifferentLabels : public std::runtime_error {
	DifferentLabels() :
		std::runtime_error("histograms have different bin limits") { }
};

Hist4D::Hist4D(
		char const* name, char const* title,
		Int_t nbins_w, Double_t const* bins_w,
		Int_t nbins_x, Double_t const* bins_x,
		Int_t nbins_y, Double_t const* bins_y,
		Int_t nbins_z, Double_t const* bins_z) :
		TNamed(name, title),
		_w_axis(new TAxis(nbins_w, bins_w)),
		_x_axis(new TAxis(nbins_x, bins_x)),
		_y_axis(new TAxis(nbins_y, bins_y)),
		_z_axis(new TAxis(nbins_z, bins_z)),
		_hists() {
	CreateAxes();
	CreateHists();
}

Hist4D::Hist4D(
		char const* name, char const* title,
		Int_t nbins_w, Double_t w_lower, Double_t w_upper,
		Int_t nbins_x, Double_t x_lower, Double_t x_upper,
		Int_t nbins_y, Double_t y_lower, Double_t y_upper,
		Int_t nbins_z, Double_t z_lower, Double_t z_upper) :
		TNamed(name, title),
		_w_axis(new TAxis(nbins_w, w_lower, w_upper)),
		_x_axis(new TAxis(nbins_x, x_lower, x_upper)),
		_y_axis(new TAxis(nbins_y, y_lower, y_upper)),
		_z_axis(new TAxis(nbins_z, z_lower, z_upper)),
		_hists() {
	CreateAxes();
	CreateHists();
}

bool Hist4D::CheckConsistency(Hist4D const* h1, Hist4D const* h2) {
	if (h1 == h2) {
		return true;
	}
	Int_t nbins_w = h1->_w_axis->GetNbins();
	Int_t nbins_x = h1->_x_axis->GetNbins();
	Int_t nbins_y = h1->_y_axis->GetNbins();
	Int_t nbins_z = h1->_z_axis->GetNbins();
	if (nbins_w != h2->_w_axis->GetNbins()
			|| nbins_x != h2->_x_axis->GetNbins()
			|| nbins_y != h2->_y_axis->GetNbins()
			|| nbins_z != h2->_z_axis->GetNbins()) {
		throw DifferentNumberOfBins();
		return false;
	}
	if (h1->_hists.size() != h2->_hists.size()) {
		throw DifferentNumberOfBins();
		return false;
	}
	bool result = true;
	result &= CheckAxisLimits(h1->_w_axis, h2->_w_axis);
	result &= CheckAxisLimits(h1->_x_axis, h2->_x_axis);
	result &= CheckAxisLimits(h1->_y_axis, h2->_y_axis);
	result &= CheckAxisLimits(h1->_z_axis, h2->_z_axis);
	result &= CheckBinLimits(h1->_w_axis, h2->_w_axis);
	result &= CheckBinLimits(h1->_x_axis, h2->_x_axis);
	result &= CheckBinLimits(h1->_y_axis, h2->_y_axis);
	result &= CheckBinLimits(h1->_z_axis, h2->_z_axis);
	result &= CheckBinLabels(h1->_w_axis, h2->_w_axis);
	result &= CheckBinLabels(h1->_x_axis, h2->_x_axis);
	result &= CheckBinLabels(h1->_y_axis, h2->_y_axis);
	result &= CheckBinLabels(h1->_z_axis, h2->_z_axis);
	return result;
}

bool Hist4D::CheckAxisLimits(TAxis const* ax1, TAxis const* ax2) {
	Double_t first_bin = ax1->GetBinWidth(1);
	Double_t last_bin = ax1->GetBinWidth(ax1->GetNbins());
	if (!TMath::AreEqualAbs(ax1->GetXmin(), ax2->GetXmin(), first_bin * 1e-10)
			|| !TMath::AreEqualAbs(ax1->GetXmax(), ax2->GetXmax(), last_bin * 1e-10)) {
		throw DifferentAxisLimits();
		return false;
	}
	return true;
}

bool Hist4D::CheckBinLimits(TAxis const* ax1, TAxis const* ax2) {
	TArrayD const* h1arr = ax1->GetXbins();
	TArrayD const* h2arr = ax2->GetXbins();
	Int_t count = h1arr->GetSize();
	if (count != 0) {
		if (h2arr->GetSize() != count) {
			throw DifferentBinLimits();
			return false;
		} else {
			for (Int_t idx = 0; idx < count; ++idx) {
				Double_t bin_width = ax1->GetBinWidth(idx);
				if (!TMath::AreEqualAbs(h1arr->GetAt(idx), h2arr->GetAt(idx), bin_width * 1e-10)) {
					throw DifferentBinLimits();
					return false;
				}
			}
		}
	}
	return true;
}

bool Hist4D::CheckBinLabels(TAxis const* ax1, TAxis const* ax2) {
	THashList* l1 = ax1->GetLabels();
	THashList* l2 = ax2->GetLabels();
	if (!l1 && !l2) {
		return true;
	} else if (!l1 || !l2) {
		throw DifferentLabels();
		return false;
	} else {
		for (Int_t idx = 1; idx <= ax1->GetNbins(); ++idx) {
			TString label1 = ax1->GetBinLabel(idx);
			TString label2 = ax2->GetBinLabel(idx);
			if (label1 != label2) {
				throw DifferentLabels();
				return false;
			}
		}
	}
	return true;
}

void Hist4D::CreateHists() {
	// Include overflow and underflow hists.
	Int_t nbins_w = _w_axis->GetNbins();
	Int_t nbins_x = _x_axis->GetNbins();
	Int_t nbins_y = _y_axis->GetNbins();
	Int_t nbins_z = _z_axis->GetNbins();
	for (Int_t idx_x = 0; idx_x < nbins_x + 2; ++idx_x) {
		for (Int_t idx_w = 0; idx_w < nbins_w + 2; ++idx_w) {
			std::string hist_name = std::string(fName.Data())
				+ "/hist" + std::to_string(idx_x) + ":" + std::to_string(idx_w);
			TArrayD const* bins_y = _y_axis->GetXbins();
			TArrayD const* bins_z = _z_axis->GetXbins();
			bool linear = (bins_y->GetSize() == 0) && (bins_z->GetSize() == 0);
			TH2D* hist;
			if (linear) {
				hist = new TH2D(
					hist_name.c_str(), "",
					nbins_y, _y_axis->GetXmin(), _y_axis->GetXmax(),
					nbins_z, _z_axis->GetXmin(), _z_axis->GetXmax());
			} else {
				hist = new TH2D(
					hist_name.c_str(), "",
					bins_y->GetSize() - 1, bins_y->GetArray(),
					bins_z->GetSize() - 1, bins_z->GetArray());
			}
			hist->GetXaxis()->SetTickLength(0.);
			hist->GetXaxis()->SetLabelOffset(999.);
			hist->GetXaxis()->SetLabelSize(0.);
			hist->GetYaxis()->SetTickLength(0.);
			hist->GetYaxis()->SetLabelOffset(999.);
			hist->GetYaxis()->SetLabelSize(0.);
			hist->SetStats(0);
			Int_t idx = (nbins_w + 2) * idx_x + idx_w;
			_hists.push_back(hist);
		}
	}
}

void Hist4D::CreateAxes() {
	std::string axes_hist_name = fName.Data();
	axes_hist_name += "/axes_hist";
	TArrayD const* bins_w = _w_axis->GetXbins();
	TArrayD const* bins_x = _x_axis->GetXbins();
	bool linear = (bins_w->GetSize() == 0) && (bins_x->GetSize() == 0);
	if (linear) {
		_axes = new TH2D(
			axes_hist_name.c_str(), "",
			_w_axis->GetNbins(), _w_axis->GetXmin(), _w_axis->GetXmax(),
			_x_axis->GetNbins(), _x_axis->GetXmin(), _x_axis->GetXmax());
	} else {
		_axes = new TH2D(
			axes_hist_name.c_str(), "",
			bins_w->GetSize() - 1, bins_w->GetArray(),
			bins_x->GetSize() - 1, bins_x->GetArray());
	}
	_axes->SetStats(0);
	_w_axis = _axes->GetXaxis();
	_x_axis = _axes->GetYaxis();
	_axes->GetXaxis()->SetTickLength(0.);
	//_axes->GetXaxis()->SetLabelOffset(999.);
	//_axes->GetXaxis()->SetLabelSize(0.);
	_axes->GetYaxis()->SetTickLength(0.);
	//_axes->GetYaxis()->SetLabelOffset(999.);
	//_axes->GetYaxis()->SetLabelSize(0.);
}

TObject* Hist4D::Clone(char const* new_name) const {
	Hist4D* result = new Hist4D();
	result->SetName(new_name);
	result->SetTitle(fTitle);
	result->_w_axis = _w_axis;
	result->_x_axis = _x_axis;
	result->_y_axis = _y_axis;
	result->_z_axis = _z_axis;
	result->_axes = static_cast<TH2D*>(_axes->Clone());
	for (TH2D* hist : _hists) {
		result->_hists.push_back(static_cast<TH2D*>(hist->Clone()));
	}
	CheckConsistency(this, result);
	return result;
}

void Hist4D::Fill(Double_t w, Double_t x, Double_t y, Double_t z) {
	Int_t idx_w = _w_axis->FindBin(w);
	Int_t idx_x = _x_axis->FindBin(x);
	Int_t nbins_w = _w_axis->GetNbins();
	Int_t nbins_x = _x_axis->GetNbins();
	Int_t idx = (nbins_w + 2) * idx_x + idx_w;
	_hists[idx]->Fill(y, z);
}

void Hist4D::Fill(Double_t w, Double_t x, Double_t y, Double_t z, Double_t weight) {
	Int_t idx_w = _w_axis->FindBin(w);
	Int_t idx_x = _x_axis->FindBin(x);
	Int_t nbins_w = _w_axis->GetNbins();
	Int_t nbins_x = _x_axis->GetNbins();
	Int_t idx = (nbins_w + 2) * idx_x + idx_w;
	_hists[idx]->Fill(y, z, weight);
}

Double_t Hist4D::GetEntries() {
	Double_t result = 0;
	for (TH2D* hist : _hists) {
		result += hist->GetEntries();
	}
	return result;
}

void Hist4D::Divide(Hist4D* other) {
	CheckConsistency(this, other);
	for (std::size_t idx = 0; idx < _hists.size(); ++idx) {
		_hists[idx]->Divide(other->_hists[idx]);
	}
}

void Hist4D::Scale(Double_t c) {
	for (TH2D* hist : _hists) {
		hist->Scale(c);
	}
}

TH2D* Hist4D::ProjectionYZ(char const* pname) {
	TArrayD const* bins_y = _y_axis->GetXbins();
	TArrayD const* bins_z = _z_axis->GetXbins();
	bool linear = (bins_y->GetSize() == 0) && (bins_z->GetSize() == 0);
	TH2D* hist_proj;
	if (linear) {
		hist_proj = new TH2D(
			pname, "",
			_y_axis->GetNbins(), _y_axis->GetXmin(), _y_axis->GetXmax(),
			_z_axis->GetNbins(), _z_axis->GetXmin(), _z_axis->GetXmax());
	} else {
		hist_proj = new TH2D(
			pname, "",
			bins_y->GetSize() - 1, bins_y->GetArray(),
			bins_z->GetSize() - 1, bins_z->GetArray());
	}
	for (TH2D* hist : _hists) {
		hist_proj->Add(hist);
	}
	return hist_proj;
}

TH2D* Hist4D::ProjectionWX(
		char const* pname,
		Int_t firstbiny, Int_t lastbiny,
		Int_t firstbinz, Int_t lastbinz) {
	TArrayD const* bins_w = _w_axis->GetXbins();
	TArrayD const* bins_x = _x_axis->GetXbins();
	Int_t nbins_w = _w_axis->GetNbins();
	Int_t nbins_x = _x_axis->GetNbins();
	bool linear = (bins_w->GetSize() == 0) && (bins_x->GetSize() == 0);
	TH2D* hist_proj;
	if (linear) {
		hist_proj = new TH2D(
			pname, "",
			_w_axis->GetNbins(), _w_axis->GetXmin(), _w_axis->GetXmax(),
			_x_axis->GetNbins(), _x_axis->GetXmin(), _x_axis->GetXmax());
	} else {
		hist_proj = new TH2D(
			pname, "",
			bins_w->GetSize() - 1, bins_w->GetArray(),
			bins_x->GetSize() - 1, bins_x->GetArray());
	}
	hist_proj->Sumw2(kFALSE);
	for (Int_t idx_x = 0; idx_x < nbins_x + 2; ++idx_x) {
		for (Int_t idx_w = 0; idx_w < nbins_w + 2; ++idx_w) {
			Int_t idx = (nbins_w + 2) * idx_x + idx_w;
			TH2D* hist = _hists[idx];
			Double_t error;
			Double_t integral = hist->IntegralAndError(
				firstbiny, lastbiny, firstbinz, lastbinz, error);
			Int_t bin = hist_proj->GetBin(idx_w, idx_x);
			hist_proj->SetBinContent(bin, integral);
			hist_proj->SetBinError(bin, error);
		}
	}
	hist_proj->Sumw2(kTRUE);
	return hist_proj;
}

void Hist4D::DrawPad(TVirtualPad* pad, char const* options) {
	// Update the axes on all histograms.
	for (TH2D* hist : _hists) {
		TArrayD const* bins_y = _y_axis->GetXbins();
		TArrayD const* bins_z = _z_axis->GetXbins();
		bool linear = (bins_y->GetSize() == 0) && (bins_z->GetSize() == 0);
		if (linear) {
			hist->SetBins(
				_y_axis->GetNbins(), _y_axis->GetXmin(), _y_axis->GetXmax(),
				_z_axis->GetNbins(), _z_axis->GetXmin(), _z_axis->GetXmax());
		} else {
			hist->SetBins(
				bins_y->GetSize() - 1, bins_y->GetArray(),
				bins_z->GetSize() - 1, bins_z->GetArray());
		}
	}
	// Find minimum and maximum.
	Double_t min = _minimum;
	Double_t max = _maximum;
	if (TMath::IsNaN(min)) {
		min = TMath::Infinity();
		for (TH2D* hist : _hists) {
			Double_t next_min = hist->GetBinContent(hist->GetMinimumBin());
			if (next_min < min) {
				min = next_min;
			}
		}
	}
	if (TMath::IsNaN(max)) {
		max = -TMath::Infinity();
		for (TH2D* hist : _hists) {
			Double_t next_max = hist->GetBinContent(hist->GetMaximumBin());
			if (next_max > max) {
				max = next_max;
			}
		}
	}
	// Draw axes.
	pad->cd();
	_axes->SetTitle(fTitle);
	_axes->SetMinimum(min);
	_axes->SetMaximum(max);
	std::string axes_options = "colz ";
	axes_options += options;
	_axes->Draw(axes_options.c_str());
	pad->Paint();
	// Create sub-pads for individual plots.
	Int_t nbins_w = _w_axis->GetNbins();
	Int_t nbins_x = _x_axis->GetNbins();
	Int_t nbins_y = _y_axis->GetNbins();
	Int_t nbins_z = _z_axis->GetNbins();
	for (Int_t idx_x = 1; idx_x < nbins_x + 1; ++idx_x) {
		for (Int_t idx_w = 1; idx_w < nbins_w + 1; ++idx_w) {
			Int_t idx = (nbins_w + 2) * idx_x + idx_w;
			std::string pad_name = pad->GetName();
			pad_name += "/subpad" + std::to_string(idx_x) + ":" + std::to_string(idx_w);
			pad->cd();
			TPad* subpad = dynamic_cast<TPad*>(
				pad->GetPrimitive(pad_name.c_str()));
			if (subpad == nullptr) {
				Double_t width = pad->GetWw() * pad->GetAbsWNDC();
				Double_t height = pad->GetWh() * pad->GetAbsHNDC();
				Double_t pad_x_lo = pad->XtoPixel(pad->XtoPad(_w_axis->GetBinLowEdge(idx_w))) / width;
				Double_t pad_y_lo = pad->YtoPixel(pad->YtoPad(_x_axis->GetBinLowEdge(idx_x))) / height;
				Double_t pad_x_hi = pad->XtoPixel(pad->XtoPad(_w_axis->GetBinUpEdge(idx_w))) / width;
				Double_t pad_y_hi = pad->YtoPixel(pad->YtoPad(_x_axis->GetBinUpEdge(idx_x))) / height;
				Double_t pad_x_min = pad->XtoPixel(pad->GetUxmin()) / width;
				Double_t pad_y_min = pad->YtoPixel(pad->GetUymin()) / height;
				Double_t pad_x_max = pad->XtoPixel(pad->GetUxmax()) / width;
				Double_t pad_y_max = pad->YtoPixel(pad->GetUymax()) / height;
				if (!(pad_x_min < pad_x_max)) {
					Double_t temp = pad_x_min;
					pad_x_min = pad_x_max;
					pad_x_max = temp;
				}
				if (!(pad_y_min < pad_y_max)) {
					Double_t temp = pad_y_min;
					pad_y_min = pad_y_max;
					pad_y_max = temp;
				}
				if (!(pad_x_lo < pad_x_hi)) {
					Double_t temp = pad_x_lo;
					pad_x_lo = pad_x_hi;
					pad_x_hi = temp;
				}
				if (!(pad_y_lo < pad_y_hi)) {
					Double_t temp = pad_y_lo;
					pad_y_lo = pad_y_hi;
					pad_y_hi = temp;
				}
				if (pad_x_lo >= pad_x_min
						&& pad_y_lo >= pad_y_min
						&& pad_x_hi <= pad_x_max
						&& pad_y_hi <= pad_y_max) {
					subpad = new TPad(
						pad_name.c_str(), "",
						pad_x_lo, 1. - pad_y_hi, pad_x_hi, 1. - pad_y_lo);
					subpad->SetLeftMargin(0.01);
					subpad->SetRightMargin(0.01);
					subpad->SetTopMargin(0.01);
					subpad->SetBottomMargin(0.01);
					subpad->SetFillStyle(4000);
					subpad->SetFrameFillStyle(4000);
					subpad->SetFrameLineColor(0);
					subpad->SetFrameLineWidth(0);
					subpad->SetBorderMode(0);
					if (pad->GetLogz()) {
						subpad->SetLogz();
					}
				}
			}
			if (subpad != nullptr) {
				pad->cd();
				subpad->Draw();
				subpad->cd();
				std::string hist_options = "cola ";
				hist_options += options;
				_hists[idx]->SetMinimum(min);
				_hists[idx]->SetMaximum(max);
				_hists[idx]->Draw(hist_options.c_str());
			}
		}
	}
	pad->cd();
}

Hist4D::~Hist4D() {
  if(_w_axis) delete _w_axis;
  if(_x_axis) delete _x_axis;
  if(_y_axis) delete _y_axis;
  if(_z_axis) delete _z_axis;
  if(_axes) delete _axes;
  for(auto e : _hists)
    if(e) delete e;
  _hists.clear();
}
