// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

#pragma once

#include <TROOT.h>
#include <TNamed.h>
#include <TPad.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFrame.h>
#include <stdexcept>

// Convenience class for plotting 4d histogram. Can't derive from TH1 because
// some of the virtual methods are only designed with up to 3 dimensions in
// mind.
class Hist4D : public TNamed {
public:
	// TODO: Store a TAxis object for each of the 4 axes, then copy it into the
	// sub-histograms as they are constructed. This also allows the _axes
	// variable to be eliminated so that it is only constructed on Draw.
	// To accomplish this, chain:
	//  * PadtoX
	//  * XtoPixel
	//  * UtoPixel (then reverse it)
	// Should look like XtoPixel(XtoPad(bin)) / UtoPixel(1).
	TAxis* _w_axis = nullptr;
	TAxis* _x_axis = nullptr;
	TAxis* _y_axis = nullptr;
	TAxis* _z_axis = nullptr;
	TH2D* _axes = nullptr;
	std::vector<TH2D*> _hists;
	Double_t _minimum = TMath::QuietNaN();
	Double_t _maximum = TMath::QuietNaN();

	void CreateHists();
	void CreateAxes();

public:
	Hist4D() { }
	Hist4D(
			char const* name, char const* title,
			Int_t nbins_w, Double_t const* bins_w,
			Int_t nbins_x, Double_t const* bins_x,
			Int_t nbins_y, Double_t const* bins_y,
			Int_t nbins_z, Double_t const* bins_z);
	Hist4D(
			char const* name, char const* title,
			Int_t nbins_w, Double_t w_lower, Double_t w_upper,
			Int_t nbins_x, Double_t x_lower, Double_t x_upper,
			Int_t nbins_y, Double_t y_lower, Double_t y_upper,
			Int_t nbins_z, Double_t z_lower, Double_t z_upper);
        ~Hist4D();

	static bool CheckConsistency(Hist4D const* h1, Hist4D const* h2);
	static bool CheckAxisLimits(TAxis const* ax1, TAxis const* ax2);
	static bool CheckBinLimits(TAxis const* ax1, TAxis const* ax2);
	static bool CheckBinLabels(TAxis const* ax1, TAxis const* ax2);

	virtual TObject* Clone(char const* new_name = "_clone") const override;

	void Fill(Double_t w, Double_t x, Double_t y, Double_t z);
	void Fill(Double_t w, Double_t x, Double_t y, Double_t z, Double_t weight);

	Double_t GetEntries();

	Bool_t IsEmpty() const;

	TAxis* GetWaxis() {
		return _w_axis;
	}
	TAxis* GetXaxis() {
		return _x_axis;
	}
	TAxis* GetYaxis() {
		return _y_axis;
	}
	TAxis* GetZaxis() {
		return _z_axis;
	}

	void SetMinimum(Double_t min = TMath::QuietNaN()) {
		_minimum = min;
	}
	void SetMaximum(Double_t max = TMath::QuietNaN()) {
		_maximum = max;
	}

	void Divide(Hist4D* other);
	void Scale(Double_t c);

	TH2D* ProjectionYZ(char const* pname = "_pyz");
	TH2D* ProjectionWX(
		char const* pname = "_pwx",
		Int_t firstbiny = 0, Int_t lastbiny = -1,
		Int_t firstbinz = 0, Int_t lastbinz = -1);
	virtual void Draw(char const* options = "col") override {
		if (!gPad) {
			gROOT->MakeDefCanvas();
		}
		DrawPad(gPad, options);
	}
	void DrawPad(TVirtualPad* pad, char const* options = "col");

	ClassDefOverride(Hist4D, 1);
};
