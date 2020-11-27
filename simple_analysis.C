#include <TLorentzVector.h>

const double J_PSI_MASS    = 3.0969;         // GeV
const double PROTON_MASS   = 0.938272081;    // GeV
const double NEUTRON_MASS  = 0.93827231;     // GeV
const double ELECTRON_MASS = 0.000510998946; // GeV
const double DEUTERON_MASS = 1.87561294257;  // GeV

const double DETECTOR_EFFICIENCY = 1.0;
const double RUN_TIME = 60. * 60. * 24. * 50.; // s
const double LUMINOSITY = 0.5*1.2e37;          // 1/(s cm^2)
const double PB_TO_CM2 = 1e-36;                // cm^2/pb

const double CROSS_SECTION_TO_RATE = LUMINOSITY * DETECTOR_EFFICIENCY * PB_TO_CM2; // 1/(s pb)
const double CROSS_SECTION_TO_COUNT = CROSS_SECTION_TO_RATE * RUN_TIME; // 1/pb
const double CROSS_SECTION_TO_HOURLY_RATE = CROSS_SECTION_TO_RATE * 60 * 60; // 1/(hr pb)

const double DEG_TO_RAD = TMath::Pi() / 180.; // 1/deg
const double RAD_TO_DEG = 180. / TMath::Pi(); // deg

void simple_analysis(char const* input_file_name, double Mll_lower, double Mll_upper) {
	TFile* file = new TFile(input_file_name);
	TTree* events = (TTree*) file->Get("tr1");
	TFile* acceptance_file = new TFile(
		"acceptance/acceptance_solid_JPsi_electron_target315_output.root");
	TH2D* acceptance_forward_angle = (TH2D*) acceptance_file->Get(
		"acceptance_ThetaP_forwardangle");
	TH2D* acceptance_large_angle = (TH2D*) acceptance_file->Get(
		"acceptance_ThetaP_largeangle");

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
	double total_rate_brem = 0.;
	double accepted_rate = 0.;
	double accepted_rate_brem = 0.;

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

		// Proton.
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

		double Mll = (p_final_lepton_p + p_final_lepton_m).M();

		if (Mll > Mll_lower && Mll < Mll_upper) {
			accepted_rate += acceptance * rate;
			accepted_rate_brem += acceptance * rate_brem;
			total_rate += rate;
			total_rate_brem += rate_brem;
		}
	}

	std::cout << "Total rate electro:    " << total_rate << " counts/hour" << std::endl;
	std::cout << "Total rate photo:      " << total_rate_brem << " counts/hour" << std::endl;
	std::cout << "Accepted rate electro: " << accepted_rate << " counts/hour" << std::endl;
	std::cout << "Accepted rate photo:   " << accepted_rate_brem << " counts/hour" << std::endl;
}

