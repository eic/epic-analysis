// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2023 Duane Byer

Double_t const PI = TMath::Pi();

Double_t const bins_x[] = {  0.01, 0.0316, 0.1, 0.316, 1. };
Double_t const bins_Q2[] = { 1., 10., 100., 1000. };
Double_t const bins_z[] = { 0.2, 0.4, 0.6, 0.8, 1. };
Double_t const bins_pt[] = { 0., 100000. };

Int_t const num_bins_x = sizeof(bins_x) / sizeof(Double_t) - 1;
Int_t const num_bins_Q2 = sizeof(bins_Q2) / sizeof(Double_t) - 1;
Int_t const num_bins_z = sizeof(bins_z) / sizeof(Double_t) - 1;
Int_t const num_bins_pt = sizeof(bins_pt) / sizeof(Double_t) - 1;
Int_t const num_bins = num_bins_x * num_bins_Q2 * num_bins_z * num_bins_pt;

Int_t get_bin(Double_t x, Double_t const* bins, Int_t num_bins) {
	Int_t idx = -1;
	for (Int_t bin = 0; bin < num_bins + 1; ++bin) {
		if (x >= bins[bin]) {
			idx += 1;
		} else {
			return idx;
		}
	}
	return idx;
}

void binning(char const* file_name, char const* output="") {
	std::cout << "Opening files." << std::endl;
	TFile* file = new TFile(file_name, "OPEN");
	TTree* events = file->Get<TTree>("tree");
	TFile* file_plots = nullptr;
	if (std::string(output) != std::string("")) {
		file_plots = new TFile(output, "RECREATE");
	}

	std::cout << "Initializing histograms." << std::endl;
	std::vector<TH2D*> angle_up_hists;
	std::vector<TH2D*> angle_down_hists;
	std::vector<TH2D*> asym_hists;
	for (Int_t bin = 0; bin < num_bins; ++bin) {
		std::string up_name = "AngleUpHist" + std::to_string(bin);
		std::string down_name = "AngleDownHist" + std::to_string(bin);
		angle_up_hists.push_back(new TH2D(
			up_name.c_str(), up_name.c_str(),
			8, -PI, PI,
			8, -PI, PI));
		angle_down_hists.push_back(new TH2D(
			down_name.c_str(), down_name.c_str(),
			8, -PI, PI,
			8, -PI, PI));
		asym_hists.push_back(nullptr);
	}

	std::cout << "Creating branch addresses." << std::endl;
	Double_t x;
	Double_t Q2;
	Double_t z;
	Double_t pt;
	Double_t phi_h;
	Double_t phi_s;
	Int_t spin_idx;
	Double_t pol;
	Double_t depol_1;
	Double_t depol_2;
	Double_t depol_3;
	Double_t depol_4;
	Double_t weight;
	events->SetBranchAddress("X", &x);
	events->SetBranchAddress("Q2", &Q2);
	events->SetBranchAddress("Z", &z);
	events->SetBranchAddress("PhPerp", &pt);
	events->SetBranchAddress("PhiH", &phi_h);
	events->SetBranchAddress("PhiS", &phi_s);
	events->SetBranchAddress("Spin_idx", &spin_idx);
	events->SetBranchAddress("Pol", &pol);
	events->SetBranchAddress("Depol1", &depol_1);
	events->SetBranchAddress("Depol2", &depol_2);
	events->SetBranchAddress("Depol3", &depol_3);
	events->SetBranchAddress("Depol4", &depol_4);
	events->SetBranchAddress("Weight", &weight);

	std::cout << "Looping over events." << std::endl;
	Int_t num_events = events->GetEntries();
	for (Int_t event = 0; event < num_events; ++event) {
		events->GetEntry(event);
		Int_t bin_x = get_bin(x, bins_x, num_bins_x);
		Int_t bin_Q2 = get_bin(Q2, bins_Q2, num_bins_Q2);
		Int_t bin_z = get_bin(z, bins_z, num_bins_z);
		Int_t bin_pt = get_bin(pt, bins_pt, num_bins_pt);
		if (!(bin_x >= 0 && bin_x < num_bins_x)
				|| !(bin_Q2 >= 0 && bin_Q2 < num_bins_Q2)
				|| !(bin_z >= 0 && bin_z < num_bins_z)
				|| !(bin_pt >= 0 && bin_pt < num_bins_pt)) {
			continue;
		}
		Int_t bin =
			bin_x + num_bins_x * (
			bin_Q2 + num_bins_Q2 * (
			bin_z + num_bins_z * (
			bin_pt)));
		Double_t S = 820.;
		Double_t y = Q2 / (S * x);
		if (y > 0.05) {
			if (spin_idx == 1) {
				angle_up_hists[bin]->Fill(phi_h, phi_s, weight);
			} else if (spin_idx == -1) {
				angle_down_hists[bin]->Fill(phi_h, phi_s, weight);
			} else {
				// Shouldn't happen.
			}
		}
	}

	std::cout << "Fitting histograms." << std::endl;
	std::vector<std::vector<Double_t> > params{ {}, {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > param_errs{ {}, {}, {}, {}, {} };
	std::vector<std::string> param_names{ "Sivers", "Collins", "sin(3φ_h - φ_s)", "sin(φ_s)", "sin(2φ_h - φ_s)" };
	TF2* fit = new TF2(
		"asymmetry",
		"[0] * ([3] * sin(x - y) + [1] * [4] * sin(x + y) + [1] * [5] * sin(3 * x - y) + [2] * [6] * sin(y) + [2] * [7] * sin(2 * x - y))",
		-PI, PI,
		-PI, PI);
	fit->FixParameter(0, pol);
	// Deal with depolarization later.
	//fit->FixParameter(1, depol_1);
	fit->FixParameter(1, 1.);
	fit->FixParameter(2, 1.);
	for (std::size_t p = 0; p < params.size(); ++p) {
		fit->SetParameter(3 + p, 0.);
		fit->SetParLimits(3 + p, -1., 1.);
	}
	for (Int_t bin = 0; bin < num_bins; ++bin) {
		Int_t bin_x = bin % num_bins_x;
		Int_t bin_Q2 = (bin / num_bins_x) % num_bins_Q2;
		Int_t bin_z = (bin / (num_bins_x * num_bins_Q2)) % num_bins_z;
		Int_t bin_pt = (bin / (num_bins_x * num_bins_Q2 * num_bins_pt)) % num_bins_pt;
		Int_t count = angle_up_hists[bin]->GetEntries();
		asym_hists[bin] = static_cast<TH2D*>(angle_up_hists[bin]->GetAsymmetry(angle_down_hists[bin]));
		std::string asym_name = "AsymmetryHist" + std::to_string(bin);
		asym_hists[bin]->SetTitle(asym_name.c_str());
		asym_hists[bin]->SetName(asym_name.c_str());
		asym_hists[bin]->Fit(fit, "IDMQ");
		for (std::size_t p = 0; p < params.size(); ++p) {
			params[p].push_back(0.);
			param_errs[p].push_back(1.);
		}
		if (count != 0) {
			for (std::size_t p = 0; p < params.size(); ++p) {
				params[p].back() = fit->GetParameter(3 + p);
				param_errs[p].back() = fit->GetParError(3 + p);
			}
		}
		std::cout << "Bin " << bin << ":" << std::endl
			<< "\tCoords:  " << bin_x << ", " << bin_Q2 << ", " << bin_z << ", " << bin_pt << std::endl
			<< "\tCount:   " << count << std::endl;
		for (std::size_t p = 0; p < params.size(); ++p) {
			std::cout << param_names[p] << ": " << params[p].back() << " ± " << param_errs[p].back() << std::endl;
		}
		if (file_plots != nullptr) {
			angle_up_hists[bin]->Write();
			angle_down_hists[bin]->Write();
			asym_hists[bin]->Write();
		}
	}
	TList* axes = new TList();
	std::vector<std::string> axis_names;
	if (num_bins_x > 1) {
		axes->Add(new TAxis(num_bins_x, bins_x));
		axis_names.push_back("x");
	}
	if (num_bins_Q2 > 1) {
		axes->Add(new TAxis(num_bins_Q2, bins_Q2));
		axis_names.push_back("Q sq.");
	}
	if (num_bins_z > 1) {
		axes->Add(new TAxis(num_bins_z, bins_z));
		axis_names.push_back("z");
	}
	if (num_bins_pt > 1) {
		axes->Add(new TAxis(num_bins_pt, bins_pt));
		axis_names.push_back("pt");
	}
	if (file_plots != nullptr) {
		file_plots->WriteObject(axes, "Axes");
		file_plots->WriteObject(&axis_names, "AxisNames");
		file_plots->WriteObject(&params, "Params");
		file_plots->WriteObject(&param_names, "ParamNames");
		file_plots->WriteObject(&param_errs, "ParamErrs");
	}
	if (file_plots != nullptr) {
		file_plots->Close();
	}
	file->Close();
}

