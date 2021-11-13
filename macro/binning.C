R__LOAD_LIBRARY(Largex)

#include <TTree.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TF2.h>

#include <fstream>
#include <vector>

#include "../src/sfset/Pavia.h"

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
	PaviaWeights* pavia = new PaviaWeights();

	Double_t const PI = TMath::Pi();

	/*
	Double_t const X1 = 1.;
	Double_t const X2 = TMath::Power(10., 1. / 3.);
	Double_t const X3 = TMath::Power(10., 2. / 3.);
	Double_t const Q1 = 1.;
	Double_t const Q2 = TMath::Power(10., 1. / 2.);
	Double_t const bins_x[] = { X1*1e-4, X2*1e-4, X3*1e-4, X1*1e-3, X2*1e-3, X3*1e-3, X1*1e-2, X2*1e-2, X3*1e-2, X1*1e-1, X2*1e-1, X3*1e-1, X1 };
	Double_t const bins_Q2[] = { Q1, Q2, Q1*1e1, Q2*1e1, Q1*1e2, Q2*1e2, Q1*1e3, Q2*1e3, Q1*1e4 };
	Double_t const bins_z[] = { 0.,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0 };
	Double_t const bins_pt[] = { 0., 100000. };
	*/
	Double_t const X1 = 1.;
	Double_t const X2 = TMath::Power(10., 1./2.);
	Double_t const bins_x[] = { X1*1e-3, X2*1e-3, X1*1e-2, X2*1e-2, X1*1e-1, X2*1e-1, X1 };
	//Double_t const bins_Q2[] = { Q1, Q2, Q1*1e1, Q2*1e1, Q1*1e2, Q2*1e2, Q1*1e3 };
	//Double_t const bins_z[] = { 0.,0.2,0.4,0.6,0.8,1.0 };
	//Double_t const bins_pt[] = { 0., 100000. };
	Double_t const bins_Q2[] = { 1., 1000000. };
	Double_t const bins_z[] = { 0., 1. };
	Double_t const bins_pt[] = { 0., 1000000. };

	Int_t const num_bins_x = sizeof(bins_x) / sizeof(Double_t) - 1;
	Int_t const num_bins_Q2 = sizeof(bins_Q2) / sizeof(Double_t) - 1;
	Int_t const num_bins_z = sizeof(bins_z) / sizeof(Double_t) - 1;
	Int_t const num_bins_pt = sizeof(bins_pt) / sizeof(Double_t) - 1;
	Int_t const num_bins = num_bins_x * num_bins_Q2 * num_bins_z * num_bins_pt;

	std::cout << "Opening files." << std::endl;
	TFile* file = new TFile(file_name, "OPEN");
	TTree* events = file->Get<TTree>("tree");
	TFile* file_plots = nullptr;
	if (std::string(output) != std::string("")) {
		file_plots = new TFile(output, "RECREATE");
	}
	Double_t xsec_total = (*file->Get<std::vector<Double_t> >("XsTotal"))[0] * 1e9; // Convert mb to fb
	Double_t weight_total = (*file->Get<std::vector<Double_t> >("WeightTotal"))[0];
	Double_t lumi = 10.; // fb^-1

	std::cout << "Initializing histograms." << std::endl;
	std::vector<TH2D*> angle_up_hists;
	std::vector<TH2D*> angle_down_hists;
	std::vector<TH2D*> asym_hists;
	std::vector<std::vector<Double_t> > asym_model{ {}, {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > asym_model_lower{ {}, {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > asym_model_upper{ {}, {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > asym_model_x{ {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > params_x{ {}, {}, {}, {} };
	std::vector<Double_t> weights_tr;
	std::vector<Double_t> weights;
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
		for (std::size_t p = 0; p < asym_model.size(); ++p) {
			asym_model[p].push_back(0.);
			asym_model_lower[p].push_back(0.);
			asym_model_upper[p].push_back(0.);
		}
		for (std::size_t ix = 0; ix < asym_model_x.size(); ++ix) {
			asym_model_x[ix].push_back(0.);
			params_x[ix].push_back(0.);
		}
		weights_tr.push_back(0.);
		weights.push_back(0.);
	}

	std::cout << "Creating branch addresses." << std::endl;
	Double_t x_tr;
	Double_t Q2_tr;
	Double_t z_tr;
	Double_t pt_tr;
	Double_t x;
	Double_t Q2;
	Double_t z;
	Double_t pt;
	Double_t phi_h;
	Double_t phi_s;
	Int_t spin_idx;
	Double_t pol_t;
	Double_t depol_1;
	Double_t depol_2;
	Double_t depol_3;
	Double_t depol_4;
	Int_t had_pid;
	Double_t weight;
	events->SetBranchAddress("TrueX", &x_tr);
	events->SetBranchAddress("TrueQSq", &Q2_tr);
	events->SetBranchAddress("TrueZ", &z_tr);
	events->SetBranchAddress("TruePhPerp", &pt_tr);
	events->SetBranchAddress("X", &x);
	events->SetBranchAddress("QSq", &Q2);
	events->SetBranchAddress("Z", &z);
	events->SetBranchAddress("PhPerp", &pt);
	events->SetBranchAddress("PhiH", &phi_h);
	events->SetBranchAddress("PhiS", &phi_s);
	events->SetBranchAddress("Spin_idx", &spin_idx);
	events->SetBranchAddress("PolT", &pol_t);
	events->SetBranchAddress("Depol1", &depol_1);
	events->SetBranchAddress("Depol2", &depol_2);
	events->SetBranchAddress("Depol3", &depol_3);
	events->SetBranchAddress("Depol4", &depol_4);
	events->SetBranchAddress("HadPID", &had_pid);
	events->SetBranchAddress("Weight", &weight);

	std::cout << "Looping over events." << std::endl;
	auto rnd = new TRandom3();
	Int_t num_events = events->GetEntries();
	for (Int_t event = 0; event < num_events; ++event) {
		if (event % 1000000 == 0) {
			std::cout << "Event " << event << std::endl;
		}
		events->GetEntry(event);
		if (had_pid != 211) {
			std::cout << "PID: " << had_pid << std::endl;
		}
		Int_t bin_x = get_bin(x, bins_x, num_bins_x);
		Int_t bin_Q2 = get_bin(Q2, bins_Q2, num_bins_Q2);
		Int_t bin_z = get_bin(z, bins_z, num_bins_z);
		Int_t bin_pt = get_bin(pt, bins_pt, num_bins_pt);
		Int_t bin_x_tr = get_bin(x_tr, bins_x, num_bins_x);
		Int_t bin_Q2_tr = get_bin(Q2_tr, bins_Q2, num_bins_Q2);
		Int_t bin_z_tr = get_bin(z_tr, bins_z, num_bins_z);
		Int_t bin_pt_tr = get_bin(pt_tr, bins_pt, num_bins_pt);
		bool bin_valid = (bin_x >= 0 && bin_x < num_bins_x)
				&& (bin_Q2 >= 0 && bin_Q2 < num_bins_Q2)
				&& (bin_z >= 0 && bin_z < num_bins_z)
				&& (bin_pt >= 0 && bin_pt < num_bins_pt);
		bool bin_valid_tr = (bin_x_tr >= 0 && bin_x_tr < num_bins_x)
				&& (bin_Q2_tr >= 0 && bin_Q2_tr < num_bins_Q2)
				&& (bin_z_tr >= 0 && bin_z_tr < num_bins_z)
				&& (bin_pt_tr >= 0 && bin_pt_tr < num_bins_pt);
		Int_t bin =
			bin_x + num_bins_x * (
			bin_Q2 + num_bins_Q2 * (
			bin_z + num_bins_z * (
			bin_pt)));
		Int_t bin_tr =
			bin_x_tr + num_bins_x * (
			bin_Q2_tr + num_bins_Q2 * (
			bin_z_tr + num_bins_z * (
			bin_pt_tr)));
		Double_t S = 4. * 10. * 100.;
		Double_t y = Q2 / (S * x);
		if (y > 0.05) {
		if (bin_valid && spin_idx == 1) {
			angle_up_hists[bin]->Fill(phi_h, phi_s, weight);
		} else if (bin_valid && spin_idx == -1) {
			angle_down_hists[bin]->Fill(phi_h, phi_s, weight);
		} else {
			// Shouldn't happen.
		}
		if (bin_valid_tr && rnd->Uniform() < 0.25) {
			Double_t sivers = pavia->Sivers(211, x_tr, z_tr, Q2_tr, pt_tr);
			Double_t sivers_lower = pavia->SiversLower(211, x_tr, z_tr, Q2_tr, pt_tr);
			Double_t sivers_upper = pavia->SiversUpper(211, x_tr, z_tr, Q2_tr, pt_tr);
			//Double_t sivers = 0.2 * x_tr;
			//Double_t sivers_lower = 0.199 * x_tr;
			//Double_t sivers_upper = 0.201 * x_tr;
			if (TMath::Finite(sivers) && TMath::Abs(sivers) < 1.) {
				asym_model[0][bin_tr] += weight * sivers;
			}
			if (TMath::Finite(sivers_lower) && TMath::Abs(sivers_lower) < 1.) {
				asym_model_lower[0][bin_tr] += weight * sivers_lower;
			}
			if (TMath::Finite(sivers_upper) && TMath::Abs(sivers_upper) < 1.) {
				asym_model_upper[0][bin_tr] += weight * sivers_upper;
			}
			asym_model_x[0][bin_tr] += weight * x_tr;
			asym_model_x[1][bin_tr] += weight * Q2_tr;
			asym_model_x[2][bin_tr] += weight * z_tr;
			asym_model_x[3][bin_tr] += weight * pt_tr;
			params_x[0][bin] += weight * x;
			params_x[1][bin] += weight * Q2;
			params_x[2][bin] += weight * z;
			params_x[3][bin] += weight * pt;
			weights_tr[bin_tr] += weight;
			weights[bin] += weight;
		}
		}
	}

	// Modify errors on histograms to account for lumi.
	for (std::size_t bin1 = 0; bin1 < num_bins; ++bin1) {
		for (Int_t bin2 = 0; bin2 < angle_up_hists[bin1]->GetNcells(); ++bin2) {
			Double_t count_up = angle_up_hists[bin1]->GetBinContent(bin2) / weight_total * lumi * xsec_total;
			Double_t error_up = TMath::Sqrt(count_up);
			angle_up_hists[bin1]->SetBinContent(bin2, count_up);
			angle_up_hists[bin1]->SetBinError(bin2, error_up);

			Double_t count_down = angle_down_hists[bin1]->GetBinContent(bin2) / weight_total * lumi * xsec_total;
			Double_t error_down = TMath::Sqrt(count_down);
			angle_down_hists[bin1]->SetBinContent(bin2, count_down);
			angle_down_hists[bin1]->SetBinError(bin2, error_down);
		}
		for (std::size_t p = 0; p < asym_model.size(); ++p) {
			asym_model[p][bin1] /= weights_tr[bin1];
			asym_model_lower[p][bin1] /= weights_tr[bin1];
			asym_model_upper[p][bin1] /= weights_tr[bin1];
		}
		weights_tr[bin1] *= lumi * xsec_total / weight_total;
		weights[bin1] *= lumi * xsec_total / weight_total;
	}

	std::ofstream fout("sivers-errors.txt");
	fout << "x\tQ2\tz\tpt\tsiv\terr\tcount\terr\tcol\terr\tpretz\terr" << std::endl;
	std::cout << "Fitting histograms." << std::endl;
	std::vector<std::vector<Double_t> > params{ {}, {}, {}, {}, {} };
	std::vector<std::vector<Double_t> > param_errs{ {}, {}, {}, {}, {} };
	std::vector<std::string> param_names{ "Sivers", "Collins", "sin(3φ_h - φ_s)", "sin(φ_s)", "sin(2φ_h - φ_s)" };
	TF2* fit = new TF2(
		"asymmetry",
		"[0] * ([3] * sin(x - y) + [1] * [4] * sin(x + y) + [1] * [5] * sin(3 * x - y) + [2] * [6] * sin(y) + [2] * [7] * sin(2 * x - y))",
		-PI, PI,
		-PI, PI);
	fit->FixParameter(0, pol_t);
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

		for (std::size_t ix = 0; ix < asym_model_x.size(); ++ix) {
			asym_model_x[ix][bin] /= weight;
			params_x[ix][bin] /= weight;
			fout << asym_model_x[ix][bin] << "\t";
		}
		fout << weights[bin] << "\t" << TMath::Sqrt(weights[bin]) << "\t";

		asym_hists[bin] = static_cast<TH2D*>(angle_up_hists[bin]->GetAsymmetry(angle_down_hists[bin]));
		std::string asym_name = "AsymmetryHist" + std::to_string(bin);
		asym_hists[bin]->SetTitle(asym_name.c_str());
		asym_hists[bin]->SetName(asym_name.c_str());
		asym_hists[bin]->Fit(fit, "IDMQ");
		for (std::size_t p = 0; p < params.size(); ++p) {
			params[p].push_back(0.);
			param_errs[p].push_back(1.);
		}
		if (asym_hists[bin]->GetEntries() != 0.) {
			for (std::size_t p = 0; p < params.size(); ++p) {
				params[p].back() = fit->GetParameter(3 + p);
				param_errs[p].back() = fit->GetParError(3 + p);
			}
		}
		std::cout << "Bin " << bin << ":" << std::endl
			<< "\tCoords:    " << bin_x << ", " << bin_Q2 << ", " << bin_z << ", " << bin_pt << std::endl
			<< "\tTotal count:  " << asym_hists[bin]->GetEntries() << std::endl
			<< "\tInjected:  " << asym_model_lower[0][bin] << " < " << asym_model[0][bin] << " < " << asym_model_upper[0][bin] << std::endl
			<< "\tExtracted: " << std::endl;
		for (std::size_t p = 0; p < params.size(); ++p) {
			std::cout << "\t\t" << param_names[p] << ": " << params[p].back() << " ± " << param_errs[p].back() << std::endl;
		}
		for (std::size_t p = 0; p < params.size(); ++p) {
			fout << params[p].back() << "\t" << param_errs[p].back() << "\t";
		}
		fout << std::endl;
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
		file_plots->WriteObject(&params_x, "ParamsX");
		file_plots->WriteObject(&asym_model, "Model");
		file_plots->WriteObject(&asym_model_lower, "ModelLower");
		file_plots->WriteObject(&asym_model_upper, "ModelUpper");
		file_plots->WriteObject(&asym_model_x, "ModelX");
		file_plots->WriteObject(&param_names, "ParamNames");
		file_plots->WriteObject(&param_errs, "ParamErrs");
	}
	if (file_plots != nullptr) {
		file_plots->Close();
	}
	file->Close();
}

