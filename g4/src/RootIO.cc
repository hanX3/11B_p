#include "RootIO.hh"
#include "OutputPath.hh"
#include "SiArray.hh"

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cstddef>
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

namespace {
void FillThreadedFileTag(char* file_name, std::size_t file_name_size)
{
  time_t t = time(nullptr);
  struct tm tt;
  localtime_r(&t, &tt);

  const G4int thread_id = G4Threading::G4GetThreadId();
  if (thread_id >= 0) {
    snprintf(file_name, file_name_size, "%d%02d%02d_%02dh%02dm%02ds_t%d", tt.tm_year + 1900, tt.tm_mon + 1, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec, thread_id);
  } else {
    snprintf(file_name, file_name_size, "%d%02d%02d_%02dh%02dm%02ds_master", tt.tm_year + 1900, tt.tm_mon + 1, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec);
  }
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO::RootIO()
{
  reaction_data.Clear();
  event_data.Clear();
  track_data.Clear();
  step_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RootIO::~RootIO()
{
  if (reaction_file) {
    delete reaction_file;
    reaction_file = nullptr;
  }
  if (event_file) {
    delete event_file;
    event_file = nullptr;
  }
  if (track_file) {
    delete track_file;
    track_file = nullptr;
  }
  if (step_file) {
    delete step_file;
    step_file = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::SetRandomSeed(ULong64_t seed)
{
  random_seed = seed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::OpenReactionFile()
{
  FillThreadedFileTag(file_name, sizeof(file_name));
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  const auto output_path = HBOutputPath::MakeOutputFilePath("reaction", file_name);

  reaction_file = new TFile(output_path.string().c_str(), "RECREATE");
  if (!reaction_file) {
    G4cout << " RootIO:: problem creating the ROOT TFile!!!" << G4endl;
    return;
  }
  G4cout << " RootIO:: successful creating the " << output_path.string() << "  !!!" << G4endl;

  run_info_tree = new TTree("RunInfo", "run information");
  run_info_tree->Branch("random_seed", &random_seed, "random_seed/l");
  run_info_tree->Fill();

  reaction_tree = new TTree("tr", "reaction simulation data");
  reaction_tree->Branch("event", &reaction_data.event, "event/L");
  reaction_tree->Branch("e_alpha1", &reaction_data.e_alpha1, "e_alpha1/D");
  reaction_tree->Branch("e_alpha2", &reaction_data.e_alpha2, "e_alpha2/D");
  reaction_tree->Branch("e_alpha3", &reaction_data.e_alpha3, "e_alpha3/D");
  reaction_tree->Branch("theta_lab_alpha1", &reaction_data.theta_lab_alpha1, "theta_lab_alpha1/D");
  reaction_tree->Branch("theta_lab_alpha2", &reaction_data.theta_lab_alpha2, "theta_lab_alpha2/D");
  reaction_tree->Branch("theta_lab_alpha3", &reaction_data.theta_lab_alpha3, "theta_lab_alpha3/D");
  reaction_tree->Branch("phi_lab_alpha1", &reaction_data.phi_lab_alpha1, "phi_lab_alpha1/D");
  reaction_tree->Branch("phi_lab_alpha2", &reaction_data.phi_lab_alpha2, "phi_lab_alpha2/D");
  reaction_tree->Branch("phi_lab_alpha3", &reaction_data.phi_lab_alpha3, "phi_lab_alpha3/D");
  reaction_tree->Branch("resonance_id", &reaction_data.resonance_id, "resonance_id/I");
  reaction_tree->Branch("branch_id", &reaction_data.branch_id, "branch_id/I");
  reaction_tree->Branch("reaction_channel", &reaction_data.reaction_channel, "reaction_channel/I");
  reaction_tree->Branch("background_mode", &reaction_data.background_mode, "background_mode/I");
  reaction_tree->Branch("gamma_resonance", &reaction_data.gamma_resonance, "gamma_resonance/I");
  reaction_tree->Branch("gamma_branch", &reaction_data.gamma_branch, "gamma_branch/I");
  reaction_tree->Branch("gamma_angular_mode", &reaction_data.gamma_angular_mode, "gamma_angular_mode/I");
  reaction_tree->Branch("n_prompt_gammas", &reaction_data.n_prompt_gammas, "n_prompt_gammas/I");
  reaction_tree->Branch("gamma1_energy", &reaction_data.gamma1_energy, "gamma1_energy/D");
  reaction_tree->Branch("gamma2_energy", &reaction_data.gamma2_energy, "gamma2_energy/D");
  reaction_tree->Branch("gamma1_theta_lab", &reaction_data.gamma1_theta_lab, "gamma1_theta_lab/D");
  reaction_tree->Branch("gamma2_theta_lab", &reaction_data.gamma2_theta_lab, "gamma2_theta_lab/D");
  reaction_tree->Branch("gamma1_phi_lab", &reaction_data.gamma1_phi_lab, "gamma1_phi_lab/D");
  reaction_tree->Branch("gamma2_phi_lab", &reaction_data.gamma2_phi_lab, "gamma2_phi_lab/D");
  reaction_tree->Branch("gamma1_theta_cm", &reaction_data.gamma1_theta_cm, "gamma1_theta_cm/D");
  reaction_tree->Branch("gamma2_theta_cm", &reaction_data.gamma2_theta_cm, "gamma2_theta_cm/D");
  reaction_tree->Branch("cos_theta_gamma_cm", &reaction_data.cos_theta_gamma_cm, "cos_theta_gamma_cm/D");
  reaction_tree->Branch("gamma_event_weight", &reaction_data.gamma_event_weight, "gamma_event_weight/D");
  reaction_tree->Branch("gamma_bias_factor", &reaction_data.gamma_bias_factor, "gamma_bias_factor/D");
  reaction_tree->Branch("event_sampling_weight", &reaction_data.event_sampling_weight, "event_sampling_weight/D");
  reaction_tree->Branch("gamma_final_state_energy_MeV", &reaction_data.gamma_final_state_energy_MeV, "gamma_final_state_energy_MeV/D");
  reaction_tree->Branch("gamma_primary_energy_MeV", &reaction_data.gamma_primary_energy_MeV, "gamma_primary_energy_MeV/D");
  reaction_tree->Branch("gamma_relative_intensity_used", &reaction_data.gamma_relative_intensity_used, "gamma_relative_intensity_used/D");
  reaction_tree->Branch("gamma_branch_fraction_used", &reaction_data.gamma_branch_fraction_used, "gamma_branch_fraction_used/D");
  reaction_tree->Branch("gamma_branch_is_upper_limit", &reaction_data.gamma_branch_is_upper_limit, "gamma_branch_is_upper_limit/I");
  reaction_tree->Branch("gamma_branch_from_relative_table", &reaction_data.gamma_branch_from_relative_table, "gamma_branch_from_relative_table/I");
  reaction_tree->Branch("gamma_cascade_generated", &reaction_data.gamma_cascade_generated, "gamma_cascade_generated/I");
  reaction_tree->Branch("enable_675_gamma_angular_distribution", &reaction_data.enable_675_gamma_angular_distribution, "enable_675_gamma_angular_distribution/I");
  reaction_tree->Branch("a1_675_gamma", &reaction_data.a1_675_gamma, "a1_675_gamma/D");
  reaction_tree->Branch("a2_675_gamma", &reaction_data.a2_675_gamma, "a2_675_gamma/D");
  reaction_tree->Branch("e_cm_p11B", &reaction_data.e_cm_p11B, "e_cm_p11B/D");
  reaction_tree->Branch("projectile_kinetic_lab", &reaction_data.projectile_kinetic_lab, "projectile_kinetic_lab/D");
  reaction_tree->Branch("projectile_px_lab", &reaction_data.projectile_px_lab, "projectile_px_lab/D");
  reaction_tree->Branch("projectile_py_lab", &reaction_data.projectile_py_lab, "projectile_py_lab/D");
  reaction_tree->Branch("projectile_pz_lab", &reaction_data.projectile_pz_lab, "projectile_pz_lab/D");
  reaction_tree->Branch("projectile_p_lab", &reaction_data.projectile_p_lab, "projectile_p_lab/D");
  reaction_tree->Branch("projectile_theta_lab", &reaction_data.projectile_theta_lab, "projectile_theta_lab/D");
  reaction_tree->Branch("projectile_phi_lab", &reaction_data.projectile_phi_lab, "projectile_phi_lab/D");
  reaction_tree->Branch("ex_max_8Be", &reaction_data.ex_max_8Be, "ex_max_8Be/D");
  reaction_tree->Branch("eaa_8Be", &reaction_data.eaa_8Be, "eaa_8Be/D");
  reaction_tree->Branch("e_alpha8Be", &reaction_data.e_alpha8Be, "e_alpha8Be/D");
  reaction_tree->Branch("cos_theta_primary_cm", &reaction_data.cos_theta_primary_cm, "cos_theta_primary_cm/D");
  reaction_tree->Branch("cos_chi_exit", &reaction_data.cos_chi_exit, "cos_chi_exit/D");
  reaction_tree->Branch("phi_primary_cm", &reaction_data.phi_primary_cm, "phi_primary_cm/D");
  reaction_tree->Branch("cos_chi_secondary_8be", &reaction_data.cos_chi_secondary_8be, "cos_chi_secondary_8be/D");
  reaction_tree->Branch("cos_theta_secondary_correlation", &reaction_data.cos_theta_secondary_correlation, "cos_theta_secondary_correlation/D");
  reaction_tree->Branch("primary_angular_mode", &reaction_data.primary_angular_mode, "primary_angular_mode/I");
  reaction_tree->Branch("enable_162_primary_angular_distribution", &reaction_data.enable_162_primary_angular_distribution, "enable_162_primary_angular_distribution/I");
  reaction_tree->Branch("a1_162_primary", &reaction_data.a1_162_primary, "a1_162_primary/D");
  reaction_tree->Branch("a2_162_primary", &reaction_data.a2_162_primary, "a2_162_primary/D");
  reaction_tree->Branch("enable_675_primary_angular_distribution", &reaction_data.enable_675_primary_angular_distribution, "enable_675_primary_angular_distribution/I");
  reaction_tree->Branch("a1_675_primary", &reaction_data.a1_675_primary, "a1_675_primary/D");
  reaction_tree->Branch("a2_675_primary", &reaction_data.a2_675_primary, "a2_675_primary/D");
  reaction_tree->Branch("primary_a1_used", &reaction_data.primary_a1_used, "primary_a1_used/D");
  reaction_tree->Branch("primary_a2_used", &reaction_data.primary_a2_used, "primary_a2_used/D");
  reaction_tree->Branch("h11b675_decay_model", &reaction_data.h11b675_decay_model, "h11b675_decay_model/I");
  reaction_tree->Branch("h11b675_decay_model_used", &reaction_data.h11b675_decay_model_used, "h11b675_decay_model_used/I");
  reaction_tree->Branch("enable_675_alpha1_secondary_angular_correlation", &reaction_data.enable_675_alpha1_secondary_angular_correlation, "enable_675_alpha1_secondary_angular_correlation/I");
  reaction_tree->Branch("h11b675_alpha_decay_model", &reaction_data.h11b675_alpha_decay_model, "h11b675_alpha_decay_model/I");
  reaction_tree->Branch("h11b675_strict_coherent_l13", &reaction_data.h11b675_strict_coherent_l13, "h11b675_strict_coherent_l13/I");
  reaction_tree->Branch("h11b675_strict_permutation_symmetrized", &reaction_data.h11b675_strict_permutation_symmetrized, "h11b675_strict_permutation_symmetrized/I");
  reaction_tree->Branch("h11b675_strict_l1_fraction", &reaction_data.h11b675_strict_l1_fraction, "h11b675_strict_l1_fraction/D");
  reaction_tree->Branch("h11b675_strict_l13_phase", &reaction_data.h11b675_strict_l13_phase, "h11b675_strict_l13_phase/D");
  reaction_tree->Branch("h11b675_strict_8be_lambda_energy_keV", &reaction_data.h11b675_strict_8be_lambda_energy_keV, "h11b675_strict_8be_lambda_energy_keV/D");
  reaction_tree->Branch("h11b675_strict_8be_reduced_width_squared_keV", &reaction_data.h11b675_strict_8be_reduced_width_squared_keV, "h11b675_strict_8be_reduced_width_squared_keV/D");
  reaction_tree->Branch("h11b675_strict_weight", &reaction_data.h11b675_strict_weight, "h11b675_strict_weight/D");
  reaction_tree->Branch("h11b675_strict_weight_max", &reaction_data.h11b675_strict_weight_max, "h11b675_strict_weight_max/D");
  reaction_tree->Branch("h11b675_strict_sampling_attempts", &reaction_data.h11b675_strict_sampling_attempts, "h11b675_strict_sampling_attempts/I");
  reaction_tree->Branch("secondary_angular_model", &reaction_data.secondary_angular_model, "secondary_angular_model/I");
  reaction_tree->Branch("secondary_a2_used", &reaction_data.secondary_a2_used, "secondary_a2_used/D");
  reaction_tree->Branch("secondary_a4_used", &reaction_data.secondary_a4_used, "secondary_a4_used/D");
  reaction_tree->Branch("background_sequential_model", &reaction_data.background_sequential_model, "background_sequential_model/I");
  reaction_tree->Branch("cos_chi_675_internal", &reaction_data.cos_chi_675_internal, "cos_chi_675_internal/D");
  reaction_tree->Branch("phi_chi_675_internal", &reaction_data.phi_chi_675_internal, "phi_chi_675_internal/D");
  reaction_tree->Branch("opening_angle_alpha12_cm", &reaction_data.opening_angle_alpha12_cm, "opening_angle_alpha12_cm/D");
  reaction_tree->Branch("opening_angle_alpha13_cm", &reaction_data.opening_angle_alpha13_cm, "opening_angle_alpha13_cm/D");
  reaction_tree->Branch("opening_angle_alpha23_cm", &reaction_data.opening_angle_alpha23_cm, "opening_angle_alpha23_cm/D");
  reaction_tree->Branch("e_alpha1_cm", &reaction_data.e_alpha1_cm, "e_alpha1_cm/D");
  reaction_tree->Branch("e_alpha2_cm", &reaction_data.e_alpha2_cm, "e_alpha2_cm/D");
  reaction_tree->Branch("e_alpha3_cm", &reaction_data.e_alpha3_cm, "e_alpha3_cm/D");
  reaction_tree->Branch("e_8be_excitation", &reaction_data.e_8be_excitation, "e_8be_excitation/D");
  reaction_tree->Branch("event_weight", &reaction_data.event_weight, "event_weight/D");
  reaction_tree->Branch("e_3alpha_cm_alpha1", &reaction_data.e_3alpha_cm_alpha1, "e_3alpha_cm_alpha1/D");
  reaction_tree->Branch("e_3alpha_cm_alpha2", &reaction_data.e_3alpha_cm_alpha2, "e_3alpha_cm_alpha2/D");
  reaction_tree->Branch("e_3alpha_cm_alpha3", &reaction_data.e_3alpha_cm_alpha3, "e_3alpha_cm_alpha3/D");
  reaction_tree->Branch("ex_8Be", &reaction_data.ex_8Be, "ex_8Be/D");
  reaction_tree->Branch("e_8Be", &reaction_data.e_8Be, "e_8Be/D");
  reaction_tree->Branch("theta_lab_8Be", &reaction_data.theta_lab_8Be, "theta_lab_8Be/D");
  reaction_tree->Branch("phi_lab_8Be", &reaction_data.phi_lab_8Be, "phi_lab_8Be/D");
  reaction_tree->Branch("x", &reaction_data.x, "x/D");
  reaction_tree->Branch("y", &reaction_data.y, "y/D");
  reaction_tree->Branch("z", &reaction_data.z, "z/D");
  reaction_tree->Branch("sigma_eval_b", &reaction_data.sigma_eval_b, "sigma_eval_b/D");
  reaction_tree->Branch("sigma_162_model_b", &reaction_data.sigma_162_model_b, "sigma_162_model_b/D");
  reaction_tree->Branch("sigma_675_model_b", &reaction_data.sigma_675_model_b, "sigma_675_model_b/D");
  reaction_tree->Branch("sigma_162_total_b", &reaction_data.sigma_162_total_b, "sigma_162_total_b/D");
  reaction_tree->Branch("sigma_675_total_b", &reaction_data.sigma_675_total_b, "sigma_675_total_b/D");
  reaction_tree->Branch("sigma_162_used_b", &reaction_data.sigma_162_used_b, "sigma_162_used_b/D");
  reaction_tree->Branch("sigma_675_used_b", &reaction_data.sigma_675_used_b, "sigma_675_used_b/D");
  reaction_tree->Branch("sigma_total_used_b", &reaction_data.sigma_total_used_b, "sigma_total_used_b/D");
  reaction_tree->Branch("sigma_162_sampling_b", &reaction_data.sigma_162_sampling_b, "sigma_162_sampling_b/D");
  reaction_tree->Branch("sigma_675_sampling_b", &reaction_data.sigma_675_sampling_b, "sigma_675_sampling_b/D");
  reaction_tree->Branch("sigma_162_directdecay_sampling_b", &reaction_data.sigma_162_directdecay_sampling_b, "sigma_162_directdecay_sampling_b/D");
  reaction_tree->Branch("sigma_675_directdecay_sampling_b", &reaction_data.sigma_675_directdecay_sampling_b, "sigma_675_directdecay_sampling_b/D");
  reaction_tree->Branch("sigma_background_sampling_b", &reaction_data.sigma_background_sampling_b, "sigma_background_sampling_b/D");
  reaction_tree->Branch("sigma_directdecay_sampling_b", &reaction_data.sigma_directdecay_sampling_b, "sigma_directdecay_sampling_b/D");
  reaction_tree->Branch("sigma_3alpha_sampling_total_b", &reaction_data.sigma_3alpha_sampling_total_b, "sigma_3alpha_sampling_total_b/D");
  reaction_tree->Branch("sigma_model_sum_b", &reaction_data.sigma_model_sum_b, "sigma_model_sum_b/D");
  reaction_tree->Branch("model_scale_factor", &reaction_data.model_scale_factor, "model_scale_factor/D");
  reaction_tree->Branch("cross_section_bias_factor", &reaction_data.cross_section_bias_factor, "cross_section_bias_factor/D");
  reaction_tree->Branch("background_bias_factor", &reaction_data.background_bias_factor, "background_bias_factor/D");
  reaction_tree->Branch("direct_decay_fraction", &reaction_data.direct_decay_fraction, "direct_decay_fraction/D");
  reaction_tree->Branch("sequential_decay_fraction_162", &reaction_data.sequential_decay_fraction_162, "sequential_decay_fraction_162/D");
  reaction_tree->Branch("sequential_decay_fraction_675", &reaction_data.sequential_decay_fraction_675, "sequential_decay_fraction_675/D");
  reaction_tree->Branch("direct_decay_fraction_162", &reaction_data.direct_decay_fraction_162, "direct_decay_fraction_162/D");
  reaction_tree->Branch("direct_decay_fraction_675", &reaction_data.direct_decay_fraction_675, "direct_decay_fraction_675/D");
  reaction_tree->Branch("enable_direct_decay", &reaction_data.enable_direct_decay, "enable_direct_decay/I");
  reaction_tree->Branch("scale_factor_162", &reaction_data.scale_factor_162, "scale_factor_162/D");
  reaction_tree->Branch("scale_factor_675", &reaction_data.scale_factor_675, "scale_factor_675/D");
  reaction_tree->Branch("sigma_background_b", &reaction_data.sigma_background_b, "sigma_background_b/D");
  reaction_tree->Branch("sigma_162_directdecay_b", &reaction_data.sigma_162_directdecay_b, "sigma_162_directdecay_b/D");
  reaction_tree->Branch("sigma_675_directdecay_b", &reaction_data.sigma_675_directdecay_b, "sigma_675_directdecay_b/D");
  reaction_tree->Branch("sigma_directdecay_b", &reaction_data.sigma_directdecay_b, "sigma_directdecay_b/D");
  reaction_tree->Branch("sigma_3alpha_eval_b", &reaction_data.sigma_3alpha_eval_b, "sigma_3alpha_eval_b/D");
  reaction_tree->Branch("sigma_gamma_162_0_b", &reaction_data.sigma_gamma_162_0_b, "sigma_gamma_162_0_b/D");
  reaction_tree->Branch("sigma_gamma_162_1_b", &reaction_data.sigma_gamma_162_1_b, "sigma_gamma_162_1_b/D");
  reaction_tree->Branch("sigma_gamma_162_total_b", &reaction_data.sigma_gamma_162_total_b, "sigma_gamma_162_total_b/D");
  reaction_tree->Branch("sigma_gamma_675_total_b", &reaction_data.sigma_gamma_675_total_b, "sigma_gamma_675_total_b/D");
  reaction_tree->Branch("sigma_gamma_total_b", &reaction_data.sigma_gamma_total_b, "sigma_gamma_total_b/D");
  reaction_tree->Branch("sigma_total_physical_all_b", &reaction_data.sigma_total_physical_all_b, "sigma_total_physical_all_b/D");
  reaction_tree->Branch("sigma_total_sampling_all_b", &reaction_data.sigma_total_sampling_all_b, "sigma_total_sampling_all_b/D");
  reaction_tree->Branch("sigma_gamma_162_0_physical_b", &reaction_data.sigma_gamma_162_0_physical_b, "sigma_gamma_162_0_physical_b/D");
  reaction_tree->Branch("sigma_gamma_162_1_physical_b", &reaction_data.sigma_gamma_162_1_physical_b, "sigma_gamma_162_1_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_physical_b", &reaction_data.sigma_gamma_675_physical_b, "sigma_gamma_675_physical_b/D");
  reaction_tree->Branch("sigma_gamma_162_0_sampling_b", &reaction_data.sigma_gamma_162_0_sampling_b, "sigma_gamma_162_0_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_162_1_sampling_b", &reaction_data.sigma_gamma_162_1_sampling_b, "sigma_gamma_162_1_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_sampling_b", &reaction_data.sigma_gamma_675_sampling_b, "sigma_gamma_675_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_ground_physical_b", &reaction_data.sigma_gamma_675_to_ground_physical_b, "sigma_gamma_675_to_ground_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_4439_physical_b", &reaction_data.sigma_gamma_675_to_4439_physical_b, "sigma_gamma_675_to_4439_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_7654_physical_b", &reaction_data.sigma_gamma_675_to_7654_physical_b, "sigma_gamma_675_to_7654_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_12710_physical_b", &reaction_data.sigma_gamma_675_to_12710_physical_b, "sigma_gamma_675_to_12710_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_15110_physical_b", &reaction_data.sigma_gamma_675_to_15110_physical_b, "sigma_gamma_675_to_15110_physical_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_ground_sampling_b", &reaction_data.sigma_gamma_675_to_ground_sampling_b, "sigma_gamma_675_to_ground_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_4439_sampling_b", &reaction_data.sigma_gamma_675_to_4439_sampling_b, "sigma_gamma_675_to_4439_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_7654_sampling_b", &reaction_data.sigma_gamma_675_to_7654_sampling_b, "sigma_gamma_675_to_7654_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_12710_sampling_b", &reaction_data.sigma_gamma_675_to_12710_sampling_b, "sigma_gamma_675_to_12710_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_675_to_15110_sampling_b", &reaction_data.sigma_gamma_675_to_15110_sampling_b, "sigma_gamma_675_to_15110_sampling_b/D");
  reaction_tree->Branch("sigma_gamma_sampling_total_b", &reaction_data.sigma_gamma_sampling_total_b, "sigma_gamma_sampling_total_b/D");
  reaction_tree->Branch("sigma_total_all_b", &reaction_data.sigma_total_all_b, "sigma_total_all_b/D");
  reaction_tree->Branch("channel_probability_162", &reaction_data.channel_probability_162, "channel_probability_162/D");
  reaction_tree->Branch("channel_probability_675", &reaction_data.channel_probability_675, "channel_probability_675/D");
  reaction_tree->Branch("channel_probability_background", &reaction_data.channel_probability_background, "channel_probability_background/D");
  reaction_tree->Branch("channel_probability_directdecay", &reaction_data.channel_probability_directdecay, "channel_probability_directdecay/D");
  reaction_tree->Branch("channel_probability_directdecay_162", &reaction_data.channel_probability_directdecay_162, "channel_probability_directdecay_162/D");
  reaction_tree->Branch("channel_probability_directdecay_675", &reaction_data.channel_probability_directdecay_675, "channel_probability_directdecay_675/D");
  reaction_tree->Branch("channel_probability_gamma", &reaction_data.channel_probability_gamma, "channel_probability_gamma/D");
  reaction_tree->Branch("probability_gamma_162_0", &reaction_data.probability_gamma_162_0, "probability_gamma_162_0/D");
  reaction_tree->Branch("probability_gamma_162_1", &reaction_data.probability_gamma_162_1, "probability_gamma_162_1/D");
  reaction_tree->Branch("probability_gamma_675_total", &reaction_data.probability_gamma_675_total, "probability_gamma_675_total/D");
  reaction_tree->Branch("probability_gamma_675_to_ground", &reaction_data.probability_gamma_675_to_ground, "probability_gamma_675_to_ground/D");
  reaction_tree->Branch("probability_gamma_675_to_4439", &reaction_data.probability_gamma_675_to_4439, "probability_gamma_675_to_4439/D");
  reaction_tree->Branch("probability_gamma_675_to_7654", &reaction_data.probability_gamma_675_to_7654, "probability_gamma_675_to_7654/D");
  reaction_tree->Branch("probability_gamma_675_to_12710", &reaction_data.probability_gamma_675_to_12710, "probability_gamma_675_to_12710/D");
  reaction_tree->Branch("probability_gamma_675_to_15110", &reaction_data.probability_gamma_675_to_15110, "probability_gamma_675_to_15110/D");
  reaction_tree->Branch("reaction", reaction_data.reaction, "reaction/C");

  if (!reaction_tree) {
    G4cout << "\n can't create tree" << G4endl;
    return;
  }
  G4cout << "\n----> Tree file is opened in " << output_path.string() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillReactionTree(H11BReactionData& data)
{
  if (!reaction_tree) return;

  reaction_data = data;
  reaction_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::CloseReactionFile()
{
  if (!reaction_file) return;

  reaction_file->cd();
  if (run_info_tree) run_info_tree->Write();
  if (reaction_tree) reaction_tree->Write();
  reaction_file->Close();
  G4cout << "\n----> reaction tree is saved.\n\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::OpenEventFile()
{
  FillThreadedFileTag(file_name, sizeof(file_name));
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  const auto output_path = HBOutputPath::MakeOutputFilePath("event", file_name);

  event_file = new TFile(output_path.string().c_str(), "RECREATE");
  if (!event_file) {
    G4cout << " RootIO:: problem creating the ROOT TFile!!!" << G4endl;
    return;
  }
  G4cout << " RootIO:: successful creating the " << output_path.string() << "  !!!" << G4endl;

  event_tree = new TTree("event", "one entry per Geant4 event");
  event_tree->Branch("event_id", &event_data.event_id, "event_id/L");

  event_tree->Branch("si_detector_id", &event_data.si_detector_id);
  event_tree->Branch("si_side", &event_data.si_side);
  event_tree->Branch("si_strip_id", &event_data.si_strip_id);
  event_tree->Branch("si_edep_MeV", &event_data.si_edep_MeV);
  event_tree->Branch("si_time_ns", &event_data.si_time_ns);

  event_tree->Branch("labr3_detector_id", &event_data.labr3_detector_id);
  event_tree->Branch("labr3_edep_MeV", &event_data.labr3_edep_MeV);
  event_tree->Branch("labr3_time_ns", &event_data.labr3_time_ns);

  event_tree->Branch("hpge_detector_id", &event_data.hpge_detector_id);
  event_tree->Branch("hpge_edep_MeV", &event_data.hpge_edep_MeV);
  event_tree->Branch("hpge_time_ns", &event_data.hpge_time_ns);

  if (!event_tree) {
    G4cout << "\n can't create event tree" << G4endl;
    return;
  }

  // Static geometry lookup used after front/back strip pairing.  Angles are
  // calculated from the nominal target centre (0,0,TargetZPos); the saved xyz
  // coordinates allow Python to recompute them for an event-specific vertex.
  si_pixel_map_tree = new TTree("si_pixel_map", "DSSD ideal-pixel centre geometry");
  si_pixel_map_tree->Branch("detector_id", &si_pixel_map_data.detector_id, "detector_id/I");
  si_pixel_map_tree->Branch("detector_model", &si_pixel_map_data.detector_model, "detector_model/I");
  si_pixel_map_tree->Branch("subarray_id", &si_pixel_map_data.subarray_id, "subarray_id/I");
  si_pixel_map_tree->Branch("module_id", &si_pixel_map_data.module_id, "module_id/I");
  si_pixel_map_tree->Branch("front_strip_id", &si_pixel_map_data.front_strip_id, "front_strip_id/I");
  si_pixel_map_tree->Branch("back_strip_id", &si_pixel_map_data.back_strip_id, "back_strip_id/I");
  si_pixel_map_tree->Branch("x_center_mm", &si_pixel_map_data.x_center_mm, "x_center_mm/D");
  si_pixel_map_tree->Branch("y_center_mm", &si_pixel_map_data.y_center_mm, "y_center_mm/D");
  si_pixel_map_tree->Branch("z_center_mm", &si_pixel_map_data.z_center_mm, "z_center_mm/D");
  si_pixel_map_tree->Branch("theta_lab_center_deg", &si_pixel_map_data.theta_lab_center_deg,
                            "theta_lab_center_deg/D");
  si_pixel_map_tree->Branch("phi_lab_center_deg", &si_pixel_map_data.phi_lab_center_deg,
                            "phi_lab_center_deg/D");

  const auto pixel_map = SiArray::BuildPixelMap();
  for (const auto& entry : pixel_map) {
    si_pixel_map_data = entry;
    si_pixel_map_tree->Fill();
  }

  G4cout << "----> event tree and " << pixel_map.size() << " Si pixel-map entries are ready in "
         << output_path.string() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillEventTree(EventData& data)
{
  if (!event_tree) return;

  event_data = data;
  event_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::CloseEventFile()
{
  if (!event_file) return;

  event_file->cd();
  if (event_tree) event_tree->Write();
  if (si_pixel_map_tree) si_pixel_map_tree->Write();
  event_file->Close();
  G4cout << "\n----> event tree and Si pixel map are saved.\n\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::OpenTrackFile()
{
  FillThreadedFileTag(file_name, sizeof(file_name));
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  const auto output_path = HBOutputPath::MakeOutputFilePath("track", file_name);

  track_file = new TFile(output_path.string().c_str(), "RECREATE");
  if (!track_file) {
    G4cout << " RootIO:: problem creating the ROOT TFile!!!" << G4endl;
    return;
  }
  G4cout << " RootIO:: successful creating the " << output_path.string() << "  !!!" << G4endl;

  track_tree = new TTree("tr", "track simulation data");
  track_tree->Branch("event", &track_data.event, "event/L");
  track_tree->Branch("track", &track_data.track, "track/I");
  track_tree->Branch("e", &track_data.e, "e/D");
  track_tree->Branch("x", &track_data.x, "x/D");
  track_tree->Branch("y", &track_data.y, "y/D");
  track_tree->Branch("z", &track_data.z, "z/D");
  track_tree->Branch("ts", &track_data.ts, "ts/D");
  track_tree->Branch("length", &track_data.length, "length/D");
  track_tree->Branch("volume", track_data.volume, "volume/C");
  track_tree->Branch("particle", track_data.particle, "particle/C");

  if (!track_tree) {
    G4cout << "\n can't create tree" << G4endl;
    return;
  }
  G4cout << "\n----> Tree file is opened in " << output_path.string() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillTrackTree(TrackData& data)
{
  if (!track_tree) return;

  track_data = data;
  track_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::CloseTrackFile()
{
  if (!track_file) return;

  track_file->cd();
  if (track_tree) track_tree->Write();
  track_file->Close();
  G4cout << "\n----> track tree is saved.\n\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::OpenStepFile()
{
  FillThreadedFileTag(file_name, sizeof(file_name));
  G4cout << "\n----> Tree file is opened in " << file_name << G4endl;

  const auto output_path = HBOutputPath::MakeOutputFilePath("step", file_name);

  step_file = new TFile(output_path.string().c_str(), "RECREATE");
  if (!step_file) {
    G4cout << " RootIO:: problem creating the ROOT TFile!!!" << G4endl;
    return;
  }
  G4cout << " RootIO:: successful creating the " << output_path.string() << "  !!!" << G4endl;

  step_tree = new TTree("tr", "step simulation data");
  step_tree->Branch("event", &step_data.event, "event/L");
  step_tree->Branch("track", &step_data.track, "track/I");
  step_tree->Branch("de", &step_data.de, "de/D");
  step_tree->Branch("pre_x", &step_data.pre_x, "pre_x/D");
  step_tree->Branch("pre_y", &step_data.pre_y, "pre_y/D");
  step_tree->Branch("pre_z", &step_data.pre_z, "pre_z/D");
  step_tree->Branch("pre_total_energy", &step_data.pre_total_energy, "pre_total_energy/D");
  step_tree->Branch("pre_kine_energy", &step_data.pre_kine_energy, "pre_kine_energy/D");
  step_tree->Branch("post_x", &step_data.post_x, "post_x/D");
  step_tree->Branch("post_y", &step_data.post_y, "post_y/D");
  step_tree->Branch("post_z", &step_data.post_z, "post_z/D");
  step_tree->Branch("post_total_energy", &step_data.post_total_energy, "post_total_energy/D");
  step_tree->Branch("post_kine_energy", &step_data.post_kine_energy, "post_kine_energy/D");
  step_tree->Branch("length", &step_data.length, "length/D");
  step_tree->Branch("volume", step_data.volume, "volume/C");
  step_tree->Branch("particle", step_data.particle, "particle/C");
  step_tree->Branch("process", step_data.process, "process/C");

  if (!step_tree) {
    G4cout << "\n can't create tree" << G4endl;
    return;
  }
  G4cout << "\n----> Tree file is opened in " << output_path.string() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::FillStepTree(StepData& data)
{
  if (!step_tree) return;

  step_data = data;
  step_tree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RootIO::CloseStepFile()
{
  if (!step_file) return;

  step_file->cd();
  if (step_tree) step_tree->Write();
  step_file->Close();
  G4cout << "\n----> step tree is saved.\n\n";
}
