// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

//R__LOAD_LIBRARY(LargexAsym)
//#include "BruAsymmetry.h"

void draw3(char const* file_name/*TString bru_dir="bruspin", TString minimizer_str="minuit"*/) {
	/*
	TFile* bin_file = new TFile(bru_dir + "/DataBinsConfig.root", "READ");
	HS::FIT::Bins* hs_bins = (HS::FIT::Bins*) bin_file->Get("HSBins");
	bin_file->Close();
	*/

	/*
	Int_t dim = hs_bins->GetNAxis();
	std::vector<TAxis> axes;
	std::vector<Int_t> num_bins;
	Int_t total_num_bins = 1;
	for (Int_t i = 0; i < dim; ++i) {
		TAxis axis = hs_bins->GetVarAxis()[i];
		axes.push_back(axis);
		num_bins.push_back(axis.GetNbins());
		total_num_bins *= num_bins.back();
	}
	*/

	/*
	std::vector<std::vector<Double_t> > params;
	std::vector<std::vector<Double_t> > param_errs;
	std::vector<TString> param_names;
	Int_t num_params = -1;
	for (Int_t b = 0; b < total_num_bins; ++b) {
		TString bin_name = hs_bins->GetBinName(i);
		TFile* result_file = new TFile(bru_dir + "/" + bin_name + "/ResultsHSMinuit2.root", "READ");

		// If the parameter names haven't been loaded yet, then load them.
		RooDataSet* param_set = result_file->Get<RooDataSet*>("FinalParameters");
		Int_t num_params_next = param_set->get()->size();
		if (num_params == -1) {
			num_params = num_params_next;
			for (Int_t p = 0; p < num_params; ++p) {
				param_names.push_back((*(param_set->get()))[i]->GetName());
			}
		}

		// Load parameter values and errors from tree.
		TTree* result_tree = result_file->Get<TTree*>("ResultTree");
		std::vector<Double_t> params_next(0., num_params);
		std::vector<Double_t> param_errs_next(0., num_params);
		for (Int_t p = 0; p < num_params; ++p) {
			result_tree->SetBranchAddress(param_names[p], &params_next[p]);
			result_tree->SetBranchAddress(param_names[p], &param_errs_next[p]);
		}
		result_tree->GetEntry(0);
		params.push_back(params_next);
		param_errs.push_back(param_errs_next);
	}
	*/

	/*
	// Transpose the parameters vectors, because they are more useful that way.
	std::vector<std::vector<Double_t> > params_new(std::vector<Double_t>(0., total_num_bins), num_params);
	std::vector<std::vector<Double_t> > param_errs_new(std::vector<Double_t>(0., total_num_bins), num_params);
	for (Int_t p = 0; p < num_params; ++p) {
		for (Int_t b = 0; b < total_num_bins; ++b) {
			params_new[p][b] = params[b][p];
			param_errs_new[p][b] = param_errs[b][p];
		}
	}
	params = std::move(params_new);
	param_errs = std::move(param_errs_new);
	*/

	TFile* file = new TFile(file_name);
	TList* axes_list = file->Get<TList>("Axes");
	auto axis_names_ptr = file->Get<std::vector<std::string> >("AxisNames");
	auto params_ptr = file->Get<std::vector<std::vector<Double_t> > >("Params");
	auto param_names_ptr = file->Get<std::vector<std::string>>("ParamNames");
	auto param_errs_ptr = file->Get<std::vector<std::vector<Double_t> > >("ParamErrs");
	auto params = *params_ptr;
	auto param_names = *param_names_ptr;
	auto param_errs = *param_errs_ptr;
	auto axis_names = *axis_names_ptr;
	std::vector<TAxis> axes;
	std::vector<Int_t> num_bins;
	Int_t total_num_bins = 1;
	Int_t num_params = param_names.size();
	for (Int_t a = 0; a < axes_list->GetEntries(); ++a) {
		TAxis axis = *static_cast<TAxis*>(axes_list->At(a));
		axes.push_back(axis);
		num_bins.push_back(axis.GetNbins());
		total_num_bins *= num_bins.back();
	}
	Int_t dim = axes.size();

	TFile* file_out = new TFile("draw3.root", "RECREATE");
	std::vector<TPad*> pads;
	for (Int_t p = 0; p < num_params; ++p) {
		// Make a plot for each parameter.
		TCanvas* pad = new TCanvas(("Canvas" + param_names[p]).c_str());
		pad->SetLogx();
		pad->SetLogy();
		pads.push_back(pad);

		if (dim == 0) {
			// For "0-d" case, no point plotting. Just output the parameter values.
			std::cout << "Parameter " << p << ": " << params[p][0] << " ± " << param_errs[p][0] << std::endl;
		} else if (dim == 1) {
			// For 1-d plots, a TGraphError will work.
			Int_t count = axes[0].GetNbins();
			std::vector<Double_t> zeroes(0., count);
			Double_t const* edges = axes[0].GetXbins()->GetArray();
			std::vector<Double_t> centers;
			for (Int_t iz = 0; iz < count; ++iz) {
				centers.push_back(0.5 * (edges[iz + 1] + edges[iz]));
			}
			Double_t const* xs = centers.data();
			Double_t const* x_errs = zeroes.data();
			Double_t const* ys = params[p].data();
			Double_t const* y_errs = param_errs[p].data();
			TGraphErrors* graph = new TGraphErrors(count, xs, ys, x_errs, y_errs);
			graph->SetTitle("");
			graph->GetXaxis()->SetTitle(axis_names[0].c_str());
			graph->GetYaxis()->SetTitle(param_names[p].c_str());
			graph->GetXaxis()->SetLimits(edges[0], edges[count]);
			graph->GetYaxis()->SetRangeUser(-0.3, 0.3);
			pad->cd();
			graph->Draw("APE");
			TLine* line_zero = new TLine(edges[0], 0., edges[count], 0.);
			line_zero->SetLineStyle(2);
			line_zero->Draw("SAME");
			if (param_names[p] == "Sivers") {
				TLine* line_target = new TLine(edges[0], 0.2, edges[count], 0.2);
				line_target->SetLineStyle(1);
				line_target->Draw("SAME");
			}
		} else if (dim == 2) {
			// Don't know what to do in the 2D case.
		} else if (dim == 3) {
			// For 3-d plots, we need a number of TGraphErrors in an array.
			TH2D* hist_axes = new TH2D(
				("Axes" + param_names[p]).c_str(), "",
				axes[0].GetNbins(), axes[0].GetXbins()->GetArray(),
				axes[1].GetNbins(), axes[1].GetXbins()->GetArray());
			hist_axes->SetStats(0);
			hist_axes->GetXaxis()->SetTitle(axis_names[0].c_str());
			hist_axes->GetYaxis()->SetTitle(axis_names[1].c_str());
			hist_axes->SetTitle(param_names[p].c_str());
			pad->cd();
			hist_axes->Draw();
			pad->Paint();
			Double_t width = pad->GetWw() * pad->GetAbsWNDC();
			Double_t height = pad->GetWh() * pad->GetAbsHNDC();
			for (Int_t iy = 1; iy < axes[1].GetNbins() + 1; ++iy) {
				for (Int_t ix = 1; ix < axes[0].GetNbins() + 1; ++ix) {
					std::string pad_name = "Pad" + std::to_string(ix) + "," + std::to_string(iy);
					Double_t x_lo = axes[0].GetBinLowEdge(ix);
					Double_t y_lo = axes[1].GetBinLowEdge(iy);
					Double_t x_hi = axes[0].GetBinUpEdge(ix);
					Double_t y_hi = axes[1].GetBinUpEdge(iy);
					Double_t pad_x_lo = pad->XtoPixel(pad->XtoPad(x_lo)) / width;
					Double_t pad_y_lo = pad->YtoPixel(pad->YtoPad(y_lo)) / height;
					Double_t pad_x_hi = pad->XtoPixel(pad->XtoPad(x_hi)) / width;
					Double_t pad_y_hi = pad->YtoPixel(pad->YtoPad(y_hi)) / height;
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
						TPad* subpad = new TPad(
							pad_name.c_str(), "",
							pad_x_lo, 1. - pad_y_hi, pad_x_hi, 1. - pad_y_lo);
						subpad->SetLeftMargin(0.);
						subpad->SetRightMargin(0.);
						subpad->SetTopMargin(0.);
						subpad->SetBottomMargin(0.);
						subpad->SetFillStyle(4000);
						subpad->SetFrameFillStyle(4000);
						pad->cd();
						subpad->Draw();
						//subpad->SetFrameLineColor(0);
						//subpad->SetFrameLineWidth(0);
						//subpad->SetBorderMode(0);

						// Draw the TGraph on the new subpad.
						Int_t count = axes[2].GetNbins();
						std::vector<Double_t> zeroes(0., count);
						Double_t const* edges = axes[2].GetXbins()->GetArray();
						std::vector<Double_t> centers;
						for (Int_t iz = 0; iz < count; ++iz) {
							centers.push_back(0.5 * (edges[iz + 1] + edges[iz]));
						}
						Double_t const* xs = centers.data();
						Double_t const* x_errs = zeroes.data();
						std::vector<Double_t> ys_vec;
						std::vector<Double_t> y_errs_vec;
						Int_t stride = axes[0].GetNbins() * axes[1].GetNbins();
						Int_t offset = (ix - 1) + axes[0].GetNbins() * (iy - 1);
						for (Int_t iz = 0; iz < count; ++iz) {
							ys_vec.push_back(params[p][offset + iz * stride]);
							y_errs_vec.push_back(param_errs[p][offset + iz * stride]);
						}
						Double_t const* ys = ys_vec.data();
						Double_t const* y_errs = y_errs_vec.data();
						TGraphErrors* graph = new TGraphErrors(count, xs, ys, x_errs, y_errs);
						graph->SetTitle("");
						graph->GetXaxis()->SetTickLength(0.);
						graph->GetXaxis()->SetLabelOffset(999.);
						graph->GetXaxis()->SetLabelSize(0.);
						graph->GetYaxis()->SetTickLength(0.);
						graph->GetYaxis()->SetLabelOffset(999.);
						graph->GetYaxis()->SetLabelSize(0.);
						graph->GetXaxis()->SetLimits(edges[0], edges[count]);
						graph->GetYaxis()->SetRangeUser(-0.3, 0.3);
						subpad->cd();
						graph->Draw("APE");
						TLine* line_zero = new TLine(edges[0], 0., edges[count], 0.);
						line_zero->SetLineStyle(2);
						line_zero->Draw("SAME");
						if (param_names[p] == "Sivers") {
							TLine* line_target = new TLine(edges[0], 0.2, edges[count], 0.2);
							line_target->SetLineStyle(1);
							line_target->Draw("SAME");
						}
					} else {
						std::cout << "Subplot off canvas." << std::endl;
					}
				}
			}
		}

		pad->Draw();
		pad->Write();
	}
	file_out->Close();
	file->Close();
}

