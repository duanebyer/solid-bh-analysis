#include <iostream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TFileInfo.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParameter.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TStyle.h>

#include "Lresolution.h"

char const* const OUTPUT_FILE_NAME = "output.root";

// FIXME
const bool DETECTOR_RESOLUTION = true;

const bool SHOW_ACCEPTANCE_PLOTS = true;
const bool SHOW_TRAJECTORY_PLOTS = true;
const bool SHOW_PHYSICS_PLOTS = true;

const double J_PSI_MASS    = 3.0969;         // GeV
const double PROTON_MASS   = 0.938272081;    // GeV
const double NEUTRON_MASS  = 0.93827231;     // GeV
const double ELECTRON_MASS = 0.000510998946; // GeV
const double DEUTERON_MASS = 1.87561294257;  // GeV

// FIXME
const double INVARIANT_MASS_MIN = 2.8; // GeV
const double INVARIANT_MASS_MAX = 3.5; // GeV

const double EG_MIN = 0.0; // GeV
const double EG_MAX = 10.0; // GeV

const double T_MIN = -1.5; // GeV
const double T_MAX = 0.0; // GeV

// Constants for converting cross-sections into rates and counts.
const double DETECTOR_EFFICIENCY = 1.0;
const double RUN_TIME = 60. * 60. * 24. * 50.; // s
const double LUMINOSITY = 0.5*1.2e37;          // 1/(s cm^2)
const double PB_TO_CM2 = 1e-36;                // cm^2/pb

const double CROSS_SECTION_TO_RATE = LUMINOSITY * DETECTOR_EFFICIENCY * PB_TO_CM2; // 1/(s pb)
const double CROSS_SECTION_TO_COUNT = CROSS_SECTION_TO_RATE * RUN_TIME; // 1/pb
const double CROSS_SECTION_TO_HOURLY_RATE = CROSS_SECTION_TO_RATE * 60 * 60; // 1/(hr pb)

const double DEG_TO_RAD = TMath::Pi() / 180.; // 1/deg
const double RAD_TO_DEG = 180. / TMath::Pi(); // deg

void draw_with_profile(TH2D* hist) {
	TCanvas* canvas = new TCanvas();
	TPad* pad_main = new TPad("main", "", 0., 0.2, 1., 1.);
	TPad* pad_profile = new TPad("profile", "", 0., 0., 1., 0.2);
	pad_main->SetLogz();
	pad_profile->SetLogy();
	TH1D* profile = hist->ProjectionX();
	pad_main->Draw();
	pad_profile->Draw();

	pad_main->cd();
	hist->SetStats(0);
	hist->Draw("colz");

	pad_profile->cd();
	profile->SetStats(0);
	profile->SetTitle("");
	profile->Draw("hist ah");
}

double sq(double x) {
	return x * x;
}

double integrate(TH1D* hist, double x_lower, double x_upper) {
	int count = 100;
	double area = 0.5 * hist->Interpolate(x_lower) + 0.5 * hist->Interpolate(x_upper);
	for (int idx = 0; idx < count; ++idx) {
		double x = x_lower + (x_upper - x_lower) / count * (idx + 0.5);
		double y = hist->Interpolate(x);
		area += y;
	}
	return area * (x_upper - x_lower) / (count + 1);
}

// TODO: Make this not take a char for the particle.
TLorentzVector scatter_by_resolution(
		TRandom* rnd,
		Lresolution resolution,
		TLorentzVector p,
		char const* particle) {
	double res_polar[4] = { 0, 0, 0, 0 };
	double p_polar[3] = { p.P(), p.Theta() * RAD_TO_DEG, p.Phi() * RAD_TO_DEG };
	int error = resolution.GetResolution(p_polar, res_polar, particle);
	if (error != 0) {
		return p;
	}
	p_polar[0] = rnd->Gaus(p_polar[0], p_polar[0] * res_polar[0]);
	p_polar[1] = rnd->Gaus(p_polar[1], res_polar[1]);
	p_polar[2] = rnd->Gaus(p_polar[2], res_polar[2]);
	TLorentzVector result = TLorentzVector();
	result.SetXYZM(
		p_polar[0] * TMath::Sin(p_polar[1] * DEG_TO_RAD) * TMath::Cos(p_polar[2] * DEG_TO_RAD),
		p_polar[0] * TMath::Sin(p_polar[1] * DEG_TO_RAD) * TMath::Sin(p_polar[2] * DEG_TO_RAD),
		p_polar[0] * TMath::Cos(p_polar[1] * DEG_TO_RAD),
		p.M());
	return result;
}

void solid_analysis(char const* input_file_name) {
	TRandom* rnd = new TRandom3();

	TFile* file = new TFile(input_file_name);
	TTree* events = (TTree*) file->Get("tr1");

	// Load the SoLID acceptance data. SoLID has two detectors: the forward
	// angle and the large angle.
	TFile* acceptance_file = new TFile(
		"acceptance/acceptance_solid_JPsi_electron_target315_output.root");
	TH2D* acceptance_forward_angle = (TH2D*) acceptance_file->Get(
		"acceptance_ThetaP_forwardangle");
	TH2D* acceptance_large_angle = (TH2D*) acceptance_file->Get(
		"acceptance_ThetaP_largeangle");
	TH2D* acceptance = (TH2D*) acceptance_forward_angle->Clone();
	acceptance->Add(acceptance_large_angle);

	// Load the detector resolution information.
	Lresolution resolution("JPsi");

	TFile* output_file = new TFile(OUTPUT_FILE_NAME, "RECREATE");

	// Plots for the final particle trajectories.
	TH2D* plot_p_proton = new TH2D(
		"p_proton", "Proton momentum",
		30, 0., 50.,
		30, 0., 6.);
	TH2D* plot_p_lepton_p = new TH2D(
		"p_lepton_p", "Positron momentum",
		30, 0., 50.,
		30, 0., 10.);
	TH2D* plot_p_lepton_m = new TH2D(
		"p_lepton_m", "Electron momentum",
		30, 0., 50.,
		30, 0., 10.);
	TH2D* plot_p_lepton_total = new TH2D(
		"p_lepton_total", "Lepton total momentum",
		30, 0., 50.,
		30, 0., 10.);
	TH2D* plot_p_initial_proton = new TH2D(
		"p_proton_initial", "Initial proton momentum",
		30, 90., 180.,
		30, 0., 2.);
	TH2D* plot_p_total = new TH2D(
		"p_total", "Total momentum",
		30, 90., 180.,
		30, 0., 4.);
	TH2D* plot_p_lepton_and_proton = new TH2D(
		"p_lepton_and_proton", "Lepton and proton momentum",
		30, 0., 10.,
		30, 0., 6.);
	plot_p_proton->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_proton->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_lepton_p->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_lepton_p->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_lepton_m->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_lepton_m->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_lepton_total->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_lepton_total->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_initial_proton->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_initial_proton->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_total->GetXaxis()->SetTitle("Angle (degrees)");
	plot_p_total->GetYaxis()->SetTitle("Momentum (GeV)");
	plot_p_lepton_and_proton->GetXaxis()->SetTitle("Lepton momentum (GeV)");
	plot_p_lepton_and_proton->GetYaxis()->SetTitle("Proton momentum (GeV)");

	TH2D* plot_pe_initial_proton = new TH2D(
		"pe_initial_proton", "Initial proton momentum/energy",
		30, 0., 2.0,
		30, 0., 1.2);
	TH2D* plot_pe_total = new TH2D(
		"pe_total", "Total momentum/energy",
		30, 0., 2.,
		30, 0., 1.);
	plot_pe_initial_proton->GetXaxis()->SetTitle("Momentum (GeV)");
	plot_pe_initial_proton->GetYaxis()->SetTitle("Energy (GeV)");
	plot_pe_total->GetXaxis()->SetTitle("Momentum (GeV)");
	plot_pe_total->GetYaxis()->SetTitle("Energy (GeV)");

	// Create the invariant mass histograms.
	TH1D* plot_Mll = new TH1D(
		"Mll", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_accepted = new TH1D(
		"Mll_accepted", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_cut = new TH1D(
		"Mll_cuts", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_brem = new TH1D(
		"Mll_brem", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_accepted_brem = new TH1D(
		"Mll_accepted_brem", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_cut_brem = new TH1D(
		"Mll_cut_brem", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_total = new TH1D(
		"Mll_total", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_accepted_total = new TH1D(
		"Mll_accepted_total", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_Mll_cut_total = new TH1D(
		"Mll_cut_total", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);

	plot_Mll->GetXaxis()->SetTitle("Invariant mass of l^{+}l^{-} (GeV)");
	plot_Mll->GetYaxis()->SetTitle("Events per hour per GeV");
	plot_Mll->SetLineColor(1);

	// FIXME
	//plot_Mll->SetMaximum(1e4);
	//plot_Mll->SetMinimum(1e-2);

	plot_Mll_total->SetLineStyle(1);
	plot_Mll_accepted_total->SetLineStyle(1);
	plot_Mll_cut_total->SetLineStyle(1);

	plot_Mll->SetLineStyle(3);
	plot_Mll_accepted->SetLineStyle(3);
	plot_Mll_cut->SetLineStyle(3);

	plot_Mll_brem->SetLineStyle(2);
	plot_Mll_accepted_brem->SetLineStyle(2);
	plot_Mll_cut_brem->SetLineStyle(2);

	plot_Mll->SetLineColor(1);
	plot_Mll_brem->SetLineColor(1);
	plot_Mll_total->SetLineColor(1);

	plot_Mll_accepted->SetLineColor(3);
	plot_Mll_accepted_brem->SetLineColor(3);
	plot_Mll_accepted_total->SetLineColor(3);

	plot_Mll_cut->SetLineColor(2);
	plot_Mll_cut_brem->SetLineColor(2);
	plot_Mll_cut_total->SetLineColor(2);

	TH1D* plot_Mpll = new TH1D(
		"Mpll", "",
		50, 4., 5.);
	plot_Mpll->GetXaxis()->SetTitle("Invariant mass of p, l^{+}, l^{-} (GeV)");
	plot_Mpll->GetYaxis()->SetTitle("Events per hour per GeV");
	// FIXME
	//plot_Mpll->SetMaximum();
	//plot_Mpll->SetMinimum();

	TH1D* plot_Eg = new TH1D(
		"Eg", "",
		50, EG_MIN, EG_MAX);
	plot_Eg->GetXaxis()->SetTitle("Photon energy E_{g} (GeV)");
	plot_Eg->GetYaxis()->SetTitle("Events per hour per GeV");
	plot_Eg->SetLineColor(1);
	// FIXME
	//plot_Eg->SetMaximum(1e2);
	//plot_Eg->SetMinimum(1e-4);

	TH1D* plot_Eg_background = new TH1D(
		"Eg_cut", "",
		50, EG_MIN, EG_MAX);
	plot_Eg_background->GetXaxis()->SetTitle("Photon energy E_{g} (GeV)");
	plot_Eg_background->GetYaxis()->SetTitle("Events per hour per GeV");
	plot_Eg_background->SetLineColor(2);

	TH1D* plot_t = new TH1D(
		"t", "",
		50, -T_MAX, -T_MIN);
	plot_t->GetXaxis()->SetTitle("-t (GeV^{2})");
	plot_t->GetYaxis()->SetTitle("Events per hour per GeV^{2}");
	// FIXME
	//plot_t->SetMaximum(1e6);
	//plot_t->SetMinimum(1e2);

	TH2D* plot_Mll_and_t = new TH2D(
		"Mll_and_t", "",
		30, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX,
		30, T_MIN, T_MAX);
	plot_Mll_and_t->GetXaxis()->SetTitle("M_{ll} (GeV)");
	plot_Mll_and_t->GetYaxis()->SetTitle("t (GeV^{2})");
	plot_Mll_and_t->GetZaxis()->SetTitle("Event rate (counts / hour)");

	TH2D* plot_Mll_and_Eg = new TH2D(
		"Mll_and_Eg", "",
		30, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX,
		30, EG_MIN, EG_MAX);
	plot_Mll_and_Eg->GetXaxis()->SetTitle("M_{ll} (GeV)");
	plot_Mll_and_Eg->GetYaxis()->SetTitle("E_{g} (GeV)");
	plot_Mll_and_Eg->GetZaxis()->SetTitle("Event rate (counts / hour)");

	TH2D* plot_Eg_and_pmom = new TH2D(
		"Eg_and_pmom", "",
		30, EG_MIN, EG_MAX,
		30, 0., 3.0);
	plot_Eg_and_pmom->GetXaxis()->SetTitle("E_{g} (GeV)");
	plot_Eg_and_pmom->GetYaxis()->SetTitle("Proton momentum (GeV)");
	plot_Eg_and_pmom->GetZaxis()->SetTitle("Event rate (counts / hour)");

	TH2D* plot_Eg_and_t = new TH2D(
		"Eg_and_t", "",
		30, EG_MIN, EG_MAX,
		30, T_MIN, T_MAX);
	plot_Eg_and_t->GetXaxis()->SetTitle("E_{g} (GeV)");
	plot_Eg_and_t->GetYaxis()->SetTitle("t (GeV^{2})");
	plot_Eg_and_t->GetZaxis()->SetTitle("Event rate (counts / hour)");

	// Create the acceptance histograms.
	TH1D* plot_acceptance_proton = new TH1D(
		"acceptance_proton", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_acceptance_lepton_p = new TH1D(
		"acceptance_lepton_p", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);
	TH1D* plot_acceptance_lepton_m = new TH1D(
		"acceptance_lepton_m", "",
		50, INVARIANT_MASS_MIN, INVARIANT_MASS_MAX);

	plot_acceptance_proton->GetXaxis()->SetTitle("Invariant mass of l^{+}l^{-} (GeV)");
	plot_acceptance_proton->GetYaxis()->SetTitle("Average acceptance");
	plot_acceptance_proton->SetLineColor(1);
	plot_acceptance_proton->SetMaximum(1e0);
	plot_acceptance_proton->SetMinimum(1e-6);

	plot_acceptance_lepton_p->GetXaxis()->SetTitle("Invariant mass of l^{+}l^{-} (GeV)");
	plot_acceptance_lepton_p->GetYaxis()->SetTitle("Average acceptance");
	plot_acceptance_lepton_p->SetLineColor(2);
	plot_acceptance_lepton_p->SetMaximum(1e0);
	plot_acceptance_lepton_p->SetMinimum(1e-6);

	plot_acceptance_lepton_m->GetXaxis()->SetTitle("Invariant mass of l^{+}l^{-} (GeV)");
	plot_acceptance_lepton_m->GetYaxis()->SetTitle("Average acceptance");
	plot_acceptance_lepton_m->SetLineColor(3);
	plot_acceptance_lepton_m->SetMaximum(1e0);
	plot_acceptance_lepton_m->SetMinimum(1e-6);

	// Plots for angles of leptons.
	TH1D* plot_lepton_theta_lab = new TH1D(
		"lepton_theta_lab", "",
		50, 0., 180.);
	TH1D* plot_lepton_proton_theta_lab = new TH1D(
		"lepton_proton_theta_lab", "",
		50, 0., 180.);
	TH1D* plot_lepton_theta_rest = new TH1D(
		"lepton_theta_rest", "",
		50, 0., 180.);
	TH1D* plot_lepton_proton_theta_rest = new TH1D(
		"lepton_proton_theta_rest", "",
		50, 0., 180.);
	TH1D* plot_lepton_theta_lab_cut = new TH1D(
		"lepton_theta_lab_cut", "",
		50, 0., 180.);
	TH1D* plot_lepton_proton_theta_lab_cut = new TH1D(
		"lepton_proton_theta_lab_cut", "",
		50, 0., 180.);
	TH1D* plot_lepton_theta_rest_cut = new TH1D(
		"lepton_theta_rest_cut", "",
		50, 0., 180.);
	TH1D* plot_lepton_proton_theta_rest_cut = new TH1D(
		"lepton_proton_theta_rest_cut", "",
		50, 0., 180.);

	plot_lepton_theta_lab->GetXaxis()->SetTitle("\\theta^{l}_{lab} (degrees)");
	plot_lepton_theta_lab->GetYaxis()->SetTitle("Event rate (counts/hour/degree)");
	plot_lepton_theta_lab->SetLineColor(1);
	plot_lepton_theta_lab_cut->SetLineColor(2);

	plot_lepton_proton_theta_lab->GetXaxis()->SetTitle("\\theta^{lp}_{lab} (degrees)");
	plot_lepton_proton_theta_lab->GetYaxis()->SetTitle("Event rate (counts/hour/degree)");
	plot_lepton_proton_theta_lab->SetLineColor(1);
	plot_lepton_proton_theta_lab_cut->SetLineColor(2);

	plot_lepton_theta_rest->GetXaxis()->SetTitle("\\theta^{l}_{rest} (degrees)");
	plot_lepton_theta_rest->GetYaxis()->SetTitle("Event rate (counts/hour/degree)");
	plot_lepton_theta_rest->SetLineColor(1);
	plot_lepton_theta_rest_cut->SetLineColor(2);

	plot_lepton_proton_theta_rest->GetXaxis()->SetTitle("\\theta^{lp}_{rest} (degrees)");
	plot_lepton_proton_theta_rest->GetYaxis()->SetTitle("Event rate (counts/hour/degree)");
	plot_lepton_proton_theta_rest->SetLineColor(1);
	plot_lepton_proton_theta_rest_cut->SetLineColor(2);

	TH2D* plot_lepton_theta_lepton_theta_lab = new TH2D(
		"lepton_theta_lepton_theta_lab", "",
		50, 0., 180.,
		50, 0., 180.);
	TH2D* plot_lepton_theta_lepton_theta_lab_cut = new TH2D(
		"lepton_theta_lepton_theta_lab", "",
		50, 0., 180.,
		50, 0., 180.);
	TH1D* plot_lepton_lepton_theta_lab = new TH1D(
		"lepton_lepton_theta_lab", "",
		50, 0., 180.);
	TH1D* plot_lepton_lepton_theta_lab_cut = new TH1D(
		"lepton_lepton_theta_lab", "",
		50, 0., 180.);
	plot_lepton_theta_lepton_theta_lab->GetXaxis()->SetTitle("l^{+} polar angle (degrees)");
	plot_lepton_theta_lepton_theta_lab->GetYaxis()->SetTitle("l^{-} polar angle (degrees)");
	plot_lepton_theta_lepton_theta_lab_cut->GetXaxis()->SetTitle("l^{+} polar angle (degrees)");
	plot_lepton_theta_lepton_theta_lab_cut->GetYaxis()->SetTitle("l^{-} polar angle (degrees)");
	plot_lepton_lepton_theta_lab->GetXaxis()->SetTitle("Angle between leptons (degrees)");
	plot_lepton_lepton_theta_lab->GetYaxis()->SetTitle("Event rate before acceptance (counts/hour)");
	plot_lepton_lepton_theta_lab_cut->GetXaxis()->SetTitle("Angle between leptons (degrees)");
	plot_lepton_lepton_theta_lab_cut->GetYaxis()->SetTitle("Event rate after acceptance (counts/hour)");

	TLorentzVector* p_initial_proton_address = nullptr;
	TLorentzVector* p_final_proton_address = nullptr;
	TLorentzVector* p_final_lepton_p_address = nullptr;
	TLorentzVector* p_final_lepton_m_address = nullptr;
	double phase_space_factor;
	double flux_factor;
	double flux_factor_brem;
	double cross_section;

	events->SetBranchAddress("L_prot_i", &p_initial_proton_address);
	events->SetBranchAddress("L_prot", &p_final_proton_address);
	events->SetBranchAddress("L_em", &p_final_lepton_m_address);
	events->SetBranchAddress("L_ep", &p_final_lepton_p_address);
	events->SetBranchAddress("psf", &phase_space_factor);
	events->SetBranchAddress("flux_factor", &flux_factor);
	events->SetBranchAddress("flux_brem", &flux_factor_brem);
	events->SetBranchAddress("crs_BH", &cross_section);

	double Eb = ((TParameter<Double_t>*) file->Get("Eb"))->GetVal();

	double total_rate = 0.;
	double total_accepted_rate = 0.;
	double total_cut_rate = 0.;

	double total_acceptance_proton = 0.;
	double total_acceptance_lepton_p = 0.;
	double total_acceptance_lepton_m = 0.;
	double total_acceptance = 0.;

	Int_t target = ((TParameter<Int_t>*) file->Get("target"))->GetVal();

	Int_t event_count = ((TParameter<Int_t>*) file->Get("Nsim"))->GetVal();
	Int_t event_count_tot = ((TParameter<Int_t>*) file->Get("Nsim_tot"))->GetVal();
	std::cout << "Event count:       " << event_count << std::endl;
	std::cout << "Event count (tot): " << event_count_tot << std::endl;
	for (Int_t i = 0; i < events->GetEntries(); ++i) {
		events->GetEntry(i);

		TLorentzVector p_initial_electron;
		p_initial_electron.SetXYZM(0., 0., Eb, ELECTRON_MASS);
		TLorentzVector p_initial_proton = *p_initial_proton_address;
		// The final state particles are filled in from the event data.
		TLorentzVector p_final_proton = *p_final_proton_address;
		TLorentzVector p_final_lepton_p = *p_final_lepton_p_address;
		TLorentzVector p_final_lepton_m = *p_final_lepton_m_address;

		double rate =
			cross_section * phase_space_factor * flux_factor *
			CROSS_SECTION_TO_HOURLY_RATE / event_count_tot;
		double rate_brem =
			cross_section * phase_space_factor * flux_factor_brem *
			CROSS_SECTION_TO_HOURLY_RATE / event_count_tot;

		// The momentum seen by the detector.
		TLorentzVector p_detected_proton = p_final_proton;
		TLorentzVector p_detected_lepton_p = p_final_lepton_p;
		TLorentzVector p_detected_lepton_m = p_final_lepton_m;

		// Get the acceptance of the detector for the final particles.
		// TODO: It might be more reasonable to use the TH2->Interpolate method.
		double acceptance_forward_angle_proton =
			acceptance_forward_angle->GetBinContent(
				acceptance_forward_angle->FindBin(
					p_final_proton.Theta() * RAD_TO_DEG,
					p_final_proton.P()));
		double acceptance_large_angle_proton =
			acceptance_large_angle->GetBinContent(
				acceptance_large_angle->FindBin(
					p_final_proton.Theta() * RAD_TO_DEG,
					p_final_proton.P()));
		// Lepton 1.
		double acceptance_forward_angle_lepton_p =
			acceptance_forward_angle->GetBinContent(
				acceptance_forward_angle->FindBin(
					p_final_lepton_p.Theta() * RAD_TO_DEG,
					p_final_lepton_p.P()));
		double acceptance_large_angle_lepton_p =
			acceptance_large_angle->GetBinContent(
				acceptance_large_angle->FindBin(
					p_final_lepton_p.Theta() * RAD_TO_DEG,
					p_final_lepton_p.P()));
		// Lepton 2.
		double acceptance_forward_angle_lepton_m =
			acceptance_forward_angle->GetBinContent(
				acceptance_forward_angle->FindBin(
					p_final_lepton_m.Theta() * RAD_TO_DEG,
					p_final_lepton_m.P()));
		double acceptance_large_angle_lepton_m =
			acceptance_large_angle->GetBinContent(
				acceptance_large_angle->FindBin(
					p_final_lepton_m.Theta() * RAD_TO_DEG,
					p_final_lepton_m.P()));

		// The acceptance is the product of all particle acceptances.
		double acceptance_proton =
			acceptance_forward_angle_proton + acceptance_large_angle_proton;
		double acceptance_lepton_p =
			acceptance_forward_angle_lepton_p + acceptance_large_angle_lepton_p;
		double acceptance_lepton_m =
			acceptance_forward_angle_lepton_m + acceptance_large_angle_lepton_m;

		double acceptance = acceptance_proton * acceptance_lepton_p * acceptance_lepton_m;

		total_acceptance_proton += acceptance_proton;
		total_acceptance_lepton_p += acceptance_lepton_p;
		total_acceptance_lepton_m += acceptance_lepton_m;
		total_acceptance += acceptance;

		// If we account for the detector resolution, then the measured
		// momentum will differ from the true momentum.
		if (DETECTOR_RESOLUTION) {
			p_detected_proton = scatter_by_resolution(
				rnd,
				resolution,
				p_final_proton,
				"p");
			p_detected_lepton_p = scatter_by_resolution(
				rnd,
				resolution,
				p_final_lepton_p,
				"e+");
			p_detected_lepton_m = scatter_by_resolution(
				rnd,
				resolution,
				p_final_lepton_m,
				"e-");
		}

		double accepted_rate = acceptance * rate;
		double accepted_rate_brem = acceptance * rate_brem;

		double Mll = (p_final_lepton_p + p_final_lepton_m).M();
		double Mll_detected = (p_detected_lepton_p + p_detected_lepton_m).M();
		double Mpll =
			(p_final_lepton_p + p_final_lepton_m + p_final_proton).M();
		double Mpll_detected =
			(p_detected_lepton_p + p_detected_lepton_m + p_detected_proton).M();
		// The 4-momentum of the virtual photon that caused the reaction.
		TLorentzVector p_photon =
			p_final_proton + p_final_lepton_p + p_final_lepton_m
			- p_initial_proton;
		TLorentzVector p_detected_final_total =
			p_detected_lepton_p + p_detected_lepton_m + p_detected_proton;
		TLorentzVector p_initial_deuteron = TLorentzVector(0., 0., 0., DEUTERON_MASS);
		TLorentzVector p_detected_diff = p_detected_final_total - p_initial_deuteron;
		double Eg = p_photon.E();
		double Eg_detected = 0.;
		if (target == 0) {
			Eg_detected =
				(p_detected_final_total - TLorentzVector(0., 0., 0., PROTON_MASS)).E();
		} else if (target == 1) {
			Eg_detected = (
				(p_detected_diff).M2() - sq(NEUTRON_MASS))
				/ (2. * (p_detected_diff.E() - p_detected_diff.Pz()));
		}
		TLorentzVector p_detected_initial_photon =
			TLorentzVector(0., 0., Eg_detected, Eg_detected);
		TLorentzVector p_total =
			p_detected_lepton_p + p_detected_lepton_m + p_detected_proton;
		TLorentzVector p_detected_initial_proton =
			p_total - p_detected_initial_photon;
		double t = (p_final_lepton_p + p_final_lepton_m - p_photon).M2();
		double t_detected =
			(p_detected_lepton_p + p_detected_lepton_m - p_detected_initial_photon).M2();

		TLorentzVector lepton_cm = p_detected_lepton_m + p_detected_lepton_p;
		TLorentzVector p_initial_electron_rest = p_initial_electron;
		TLorentzVector p_detected_lepton_p_rest = p_detected_lepton_p;
		TLorentzVector p_detected_lepton_m_rest = p_detected_lepton_m;
		TLorentzVector p_detected_proton_rest = p_detected_proton;
		p_initial_electron_rest.Boost(-lepton_cm.BoostVector());
		p_detected_lepton_p_rest.Boost(-lepton_cm.BoostVector());
		p_detected_lepton_m_rest.Boost(-lepton_cm.BoostVector());
		p_detected_proton_rest.Boost(-lepton_cm.BoostVector());

		double lepton_theta_lab = p_detected_lepton_p.Theta();
		double lepton_proton_theta_lab = TMath::ACos(p_detected_lepton_p.Vect().Unit().Dot(p_detected_proton.Vect().Unit()));
		double lepton_theta_rest = TMath::ACos(p_detected_lepton_p_rest.Vect().Unit().Dot(p_initial_electron_rest.Vect().Unit()));
		double lepton_proton_theta_rest = TMath::ACos(p_detected_lepton_p_rest.Vect().Unit().Dot(p_detected_proton_rest.Vect().Unit()));

		// FIXME
		/*if (lepton_proton_theta_rest < 40. * DEG_TO_RAD || lepton_proton_theta_rest > 140. * DEG_TO_RAD) {
			continue;
		}*/

		// Apply any cuts.
		double cut = 1.;
		double cut_brem = 1.;
		// Require that the virtual photon be "quasi-real" (nearly on-
		// shell). This should already be guaranteed by the generator.
		/*if (p_photon.M2() > 0.05) {
			// The brem case should always have p_photon.M2 = 0, so ignore it.
			cut = 0.;
		}*/
		// Require that at least one lepton can use CC. TODO: Figure out
		// what CC is.
		/*if (p_detected_lepton_p.P() > 4.9 && p_detected_lepton_m.P() > 4.9) {
			cut = 0.;
			cut_brem = 0.;
		}*/
		// Do time-of-flight cuts on the proton. For smaller momenta, the
		// proton cannot be distinguished from a kaon.
		if (acceptance_forward_angle_proton != 0. && p_detected_proton.P() > 4.4) {
			cut = 0.;
			cut_brem = 0.;
		}
		if (acceptance_large_angle_proton != 0. && p_detected_proton.P() > 2.) {
			cut = 0.;
			cut_brem = 0.;
		}
		double cut_rate = cut * accepted_rate;
		double cut_rate_brem = cut_brem * accepted_rate_brem;

		total_rate += rate;
		total_accepted_rate += accepted_rate;
		total_cut_rate += cut_rate;

		if (SHOW_ACCEPTANCE_PLOTS) {
			double bin_width = plot_Mll->GetBinWidth(0);
			plot_acceptance_proton->Fill(Mll, acceptance_proton * rate / bin_width);
			plot_acceptance_lepton_p->Fill(Mll, acceptance_lepton_p * rate / bin_width);
			plot_acceptance_lepton_m->Fill(Mll, acceptance_lepton_m * rate / bin_width);
		}
		if (SHOW_TRAJECTORY_PLOTS) {
			if (rate != 0.) {
				plot_p_proton->Fill(
					p_final_proton.Theta() * RAD_TO_DEG,
					p_final_proton.P(),
					rate);
				plot_p_lepton_p->Fill(
					p_final_lepton_p.Theta() * RAD_TO_DEG,
					p_final_lepton_p.P(),
					rate);
				plot_p_lepton_m->Fill(
					p_final_lepton_m.Theta() * RAD_TO_DEG,
					p_final_lepton_m.P(),
					rate);
				plot_p_lepton_total->Fill(
					(p_detected_lepton_m + p_detected_lepton_p).Theta() * RAD_TO_DEG,
					(p_detected_lepton_m + p_detected_lepton_p).P(),
					rate);
				plot_p_initial_proton->Fill(
					p_detected_initial_proton.Theta() * RAD_TO_DEG,
					p_detected_initial_proton.P(),
					rate);
				plot_p_total->Fill(
					(p_total - p_initial_electron).Theta() * RAD_TO_DEG,
					(p_total - p_initial_electron).P(),
					rate);
				plot_p_lepton_and_proton->Fill(
					p_detected_lepton_m.P(),
					p_detected_proton.P(),
					rate);

				plot_pe_initial_proton->Fill(
					p_detected_initial_proton.P(),
					p_detected_initial_proton.E(),
					rate);
				plot_pe_total->Fill(
					(p_total - p_initial_electron).P(),
					(p_total - p_initial_electron).E(),
					rate);
			}
		}
		if (SHOW_PHYSICS_PLOTS) {
			if (rate != 0.) {
				double bin_width_Mll = plot_Mll->GetBinWidth(0);
				plot_Mll->Fill(Mll, rate / bin_width_Mll);
				plot_Mll_brem->Fill(Mll, rate_brem / bin_width_Mll);
				plot_Mll_total->Fill(Mll, (rate + rate_brem) / bin_width_Mll);

				double bin_width_theta = plot_lepton_theta_lab->GetBinWidth(0);
				plot_lepton_theta_lab->Fill(lepton_theta_lab * RAD_TO_DEG, rate / bin_width_theta);
				plot_lepton_proton_theta_lab->Fill(lepton_proton_theta_lab * RAD_TO_DEG, rate / bin_width_theta);
				plot_lepton_theta_rest->Fill(lepton_theta_rest * RAD_TO_DEG, rate / bin_width_theta);
				plot_lepton_proton_theta_rest->Fill(lepton_proton_theta_rest * RAD_TO_DEG, rate / bin_width_theta);

				plot_lepton_lepton_theta_lab->Fill(
					TMath::ACos(p_final_lepton_m.Vect().Unit().Dot(p_final_lepton_p.Vect().Unit())) * RAD_TO_DEG,
					rate);
				plot_lepton_theta_lepton_theta_lab->Fill(p_final_lepton_p.Theta() * RAD_TO_DEG, p_final_lepton_m.Theta() * RAD_TO_DEG, rate);

				double bin_width_Eg = plot_Eg->GetBinWidth(0);
				double bin_width_t = plot_t->GetBinWidth(0);
				plot_Eg->Fill(Eg, rate / bin_width_Eg);
				plot_t->Fill(-t, rate / bin_width_t);
			}
			if (accepted_rate != 0.) {
				double bin_width = plot_Mll_accepted->GetBinWidth(0);
				plot_Mll_accepted->Fill(Mll_detected, accepted_rate / bin_width);
				plot_Mll_accepted_brem->Fill(Mll_detected, accepted_rate_brem / bin_width);
				plot_Mll_accepted_total->Fill(Mll_detected, (accepted_rate + accepted_rate_brem) / bin_width);
			}
			if (cut_rate != 0.) {
				double bin_width = plot_Mll_cut->GetBinWidth(0);
				plot_Mll_cut->Fill(Mll_detected, cut_rate / bin_width);
				plot_Mll_cut_brem->Fill(Mll_detected, cut_rate_brem / bin_width);
				plot_Mll_cut_total->Fill(Mll_detected, (cut_rate + cut_rate_brem) / bin_width);
				plot_Mll_and_Eg->Fill(Mll_detected, Eg_detected, cut_rate);
				plot_Eg_and_pmom->Fill(Eg_detected, p_detected_proton.P(), cut_rate);
				plot_Mll_and_t->Fill(Mll_detected, t_detected, cut_rate);
				plot_Eg_and_t->Fill(Eg_detected, t_detected, cut_rate);

				// TODO: Should these show after or before acceptance? (detected?)
				double bin_width_Mpll = plot_Mpll->GetBinWidth(0);
				plot_Mpll->Fill(Mpll_detected, cut_rate / bin_width_Mpll);

				double bin_width_theta = plot_lepton_theta_lab_cut->GetBinWidth(0);
				plot_lepton_theta_lab_cut->Fill(lepton_theta_lab * RAD_TO_DEG, cut_rate / bin_width_theta);
				plot_lepton_proton_theta_lab_cut->Fill(lepton_proton_theta_lab * RAD_TO_DEG, cut_rate / bin_width_theta);
				plot_lepton_theta_rest_cut->Fill(lepton_theta_rest * RAD_TO_DEG, cut_rate / bin_width_theta);
				plot_lepton_proton_theta_rest_cut->Fill(lepton_proton_theta_rest * RAD_TO_DEG, cut_rate / bin_width_theta);

				plot_lepton_lepton_theta_lab_cut->Fill(
					TMath::ACos(p_detected_lepton_m.Vect().Unit().Dot(p_detected_lepton_p.Vect().Unit())) * RAD_TO_DEG,
					cut_rate);
				plot_lepton_theta_lepton_theta_lab_cut->Fill(p_detected_lepton_p.Theta() * RAD_TO_DEG, p_detected_lepton_m.Theta() * RAD_TO_DEG, rate);

				if (Mll_detected > J_PSI_MASS - 0.06 && Mll_detected < J_PSI_MASS + 0.06) {
					double bin_width_Eg = plot_Eg_background->GetBinWidth(0);
					plot_Eg_background->Fill(Eg_detected, cut_rate / bin_width_Eg);
				}
			}
		}
	}

	plot_acceptance_proton->Divide(plot_Mll);
	plot_acceptance_lepton_p->Divide(plot_Mll);
	plot_acceptance_lepton_m->Divide(plot_Mll);

	TLine* jpsi_line = new TLine(J_PSI_MASS, 0., J_PSI_MASS, 1000.);
	if (SHOW_ACCEPTANCE_PLOTS) {
		draw_with_profile(acceptance);
		TCanvas* canvas_acceptance = new TCanvas();
		canvas_acceptance->SetLogy();
		plot_acceptance_proton->Draw("hist");
		plot_acceptance_lepton_p->Draw("hist same");
		plot_acceptance_lepton_m->Draw("hist same");
		jpsi_line->Draw();
		TLegend* legend = new TLegend(0.6, 0.1, 0.9, 0.4);
		legend->AddEntry(plot_acceptance_proton, "p acceptance", "l");
		legend->AddEntry(plot_acceptance_lepton_p, "l^{+} acceptance", "l");
		legend->AddEntry(plot_acceptance_lepton_m, "l^{-} acceptance", "l");
		legend->Draw();

		plot_acceptance_proton->Write();
		plot_acceptance_lepton_p->Write();
		plot_acceptance_lepton_m->Write();
	}
	if (SHOW_TRAJECTORY_PLOTS) {
		draw_with_profile(plot_p_proton);
		draw_with_profile(plot_p_lepton_p);
		draw_with_profile(plot_p_lepton_m);
		draw_with_profile(plot_p_lepton_total);
		draw_with_profile(plot_p_initial_proton);
		draw_with_profile(plot_p_total);

		TCanvas* canvas_p_lepton_and_proton = new TCanvas();
		canvas_p_lepton_and_proton->SetLogz();
		plot_p_lepton_and_proton->Draw("colz");

		TCanvas* canvas_pe_initial_proton = new TCanvas();
		canvas_pe_initial_proton->SetLogz();
		plot_pe_initial_proton->Draw("colz");

		TCanvas* canvas_pe_total = new TCanvas();
		canvas_pe_total->SetLogz();
		plot_pe_total->Draw("colz");

		plot_p_proton->Write();
		plot_p_lepton_p->Write();
		plot_p_lepton_m->Write();
		plot_p_lepton_total->Write();
		plot_p_initial_proton->Write();
		plot_p_total->Write();
		plot_p_lepton_and_proton->Write();
		plot_pe_initial_proton->Write();
		plot_pe_total->Write();
	}
	if (SHOW_PHYSICS_PLOTS) {
		// Show Mll plots.
		TCanvas* canvas_Mll = new TCanvas();
		canvas_Mll->SetLogy();
		plot_Mll->Draw("hist");
		plot_Mll_brem->Draw("hist same");
		plot_Mll_total->Draw("hist same");
		plot_Mll_accepted->Draw("hist same");
		plot_Mll_accepted_brem->Draw("hist same");
		plot_Mll_accepted_total->Draw("hist same");
		plot_Mll_cut->Draw("hist same");
		plot_Mll_cut_brem->Draw("hist same");
		plot_Mll_cut_total->Draw("hist same");
		jpsi_line->Draw();

		TCanvas* canvas_Mpll = new TCanvas();
		canvas_Mpll->SetLogy();
		plot_Mpll->Draw("hist");

		TCanvas* canvas_Eg = new TCanvas();
		canvas_Eg->SetLogy();
		plot_Eg->Draw("hist");
		plot_Eg_background->Draw("hist same");

		TCanvas* canvas_t = new TCanvas();
		canvas_t->SetLogy();
		plot_t->Draw("hist");

		TCanvas* canvas_Mll_and_t = new TCanvas();
		canvas_Mll_and_t->SetLogz();
		plot_Mll_and_t->Draw("colz");

		TCanvas* canvas_Eg_and_t = new TCanvas();
		canvas_Eg_and_t->SetLogz();
		plot_Eg_and_t->Draw("colz");

		TCanvas* canvas_Mll_and_Eg = new TCanvas();
		canvas_Mll_and_Eg->SetLogz();
		plot_Mll_and_Eg->Draw("colz");

		TCanvas* canvas_Eg_and_pmom = new TCanvas();
		canvas_Eg_and_pmom->SetLogz();
		plot_Eg_and_pmom->Draw("colz");

		TCanvas* canvas_lepton_theta_lab = new TCanvas();
		plot_lepton_theta_lab->Draw("hist");
		plot_lepton_theta_lab_cut->Draw("hist same");

		TCanvas* canvas_lepton_proton_theta_lab = new TCanvas();
		plot_lepton_proton_theta_lab->Draw("hist");
		plot_lepton_proton_theta_lab_cut->Draw("hist same");

		TCanvas* canvas_lepton_theta_rest = new TCanvas();
		plot_lepton_theta_rest->Draw("hist");
		plot_lepton_theta_rest_cut->Draw("hist same");

		TCanvas* canvas_lepton_proton_theta_rest = new TCanvas();
		plot_lepton_proton_theta_rest->Draw("hist");
		plot_lepton_proton_theta_rest_cut->Draw("hist same");

		TCanvas* canvas_lepton_lepton_theta_lab = new TCanvas();
		plot_lepton_lepton_theta_lab->Draw("hist");
		TCanvas* canvas_lepton_lepton_theta_lab_cut = new TCanvas();
		plot_lepton_lepton_theta_lab_cut->Draw("hist");

		TCanvas* canvas_lepton_theta_lepton_theta_lab = new TCanvas();
		plot_lepton_theta_lepton_theta_lab->Draw("colz");
		TCanvas* canvas_lepton_theta_lepton_theta_lab_cut = new TCanvas();
		plot_lepton_theta_lepton_theta_lab_cut->Draw("colz");

		plot_Mll->Write();
		plot_Mll_brem->Write();
		plot_Mll_accepted->Write();
		plot_Mll_accepted_brem->Write();
		plot_Mll_cut->Write();
		plot_Mll_cut_brem->Write();
		plot_Mpll->Write();
		plot_Eg->Write();
		plot_Eg_background->Write();
		plot_t->Write();

		plot_Mll_and_t->Write();
		plot_Mll_and_Eg->Write();
		plot_Eg_and_pmom->Write();
		plot_Eg_and_t->Write();
	}

	std::cout << "Total rate:          " << total_rate << std::endl;
	std::cout << "Total accepted rate: " << total_accepted_rate << std::endl;
	std::cout << "Total cut rate:      " << total_cut_rate << std::endl;
	std::cout << "Fraction accepted:   " << total_cut_rate / total_rate << std::endl;
	std::cout << std::endl;

	std::cout << "Average proton acceptance:          "
		<< (100 * total_acceptance_proton / event_count_tot) << '%' << std::endl;
	std::cout << "Average lepton+ acceptance:         "
		<< (100 * total_acceptance_lepton_p / event_count_tot) << '%' << std::endl;
	std::cout << "Average lepton+ acceptance:         "
		<< (100 * total_acceptance_lepton_m / event_count_tot) << '%' << std::endl;
	std::cout << "Total cross section: " << plot_Mll->Integral() << std::endl;

	double lower = J_PSI_MASS - 0.06;
	double upper = J_PSI_MASS + 0.06;
	std::cout << "Integrated rates about J/psi peak by 120 MeV" << std::endl;
	std::cout << "Before acceptance:      " << integrate(plot_Mll, lower, upper) << " /hour" << std::endl;
	std::cout << "Before acceptance brem: " << integrate(plot_Mll_brem, lower, upper) << " /hour" << std::endl;
	std::cout << "Background rate:        " << integrate(plot_Mll_cut, lower, upper) << " /hour" << std::endl;
	std::cout << "Background rate brem:   " << integrate(plot_Mll_cut_brem, lower, upper) << " /hour" << std::endl;

	for (int i = 0; i < 13; ++i) {
		lower = i * 0.2 + EG_MIN;
		upper = lower + 0.2;
		std::cout << "Background rate from " << lower << " to " << upper << ": " << integrate(plot_Eg_background, lower, upper) << std::endl;
	}
	std::cout << "Below threshold total rate: " << integrate(plot_Eg_background, EG_MIN, 8.2) << std::endl;
	std::cout << "Above threshold total rate: " << integrate(plot_Eg_background, 8.2, EG_MAX) << std::endl;

	std::cout << "Eg plot integral: " << integrate(plot_Eg, EG_MIN, EG_MAX) << std::endl;
}

