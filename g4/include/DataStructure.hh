#ifndef DataStructure_H
#define DataStructure_H 1

#include <globals.hh>
#include <cstring>
#include <limits>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct H11BReactionData
{
  G4long event;
  G4double e_alpha1;
  G4double e_alpha2;
  G4double e_alpha3;
  G4double theta_lab_alpha1;
  G4double theta_lab_alpha2;
  G4double theta_lab_alpha3;
  G4double phi_lab_alpha1;
  G4double phi_lab_alpha2;
  G4double phi_lab_alpha3;

  G4int resonance_id; // 162 or 675
  G4int branch_id;    // 0: alpha0, 1: alpha1

  G4int reaction_channel; // 0: Resonance162, 1: Resonance675, 2: DirectDecay3Alpha, 3: GammaCapture12C
  G4int background_mode;  // 2: PhaseSpace; kept for directdecay backward compatibility

  G4int gamma_resonance;    // 0: none, 162, 675
  G4int gamma_branch;       // 0: none, 1: 162 gamma0, 2: 162 gamma1 cascade, 3: 675 effective
  G4int gamma_angular_mode; // 0: none, 1: isotropic, 2: fixed A1/A2
  G4int n_prompt_gammas;
  G4double gamma1_energy;
  G4double gamma2_energy;
  G4double gamma1_theta_lab;
  G4double gamma2_theta_lab;
  G4double gamma1_phi_lab;
  G4double gamma2_phi_lab;
  G4double gamma1_theta_cm;
  G4double gamma2_theta_cm;
  G4double cos_theta_gamma_cm;
  G4double gamma_event_weight;
  G4double gamma_bias_factor;
  G4double event_sampling_weight;
  G4double gamma_final_state_energy_MeV;
  G4double gamma_primary_energy_MeV;
  G4double gamma_relative_intensity_used;
  G4double gamma_branch_fraction_used;
  G4int gamma_branch_is_upper_limit;
  G4int gamma_branch_from_relative_table;
  G4int gamma_cascade_generated;
  G4int enable_675_gamma_angular_distribution;
  G4double a1_675_gamma;
  G4double a2_675_gamma;

  G4double e_cm_p11B;
  G4double ex_max_8Be;
  G4double eaa_8Be;
  G4double e_alpha8Be;

  G4double cos_theta_primary_cm;
  G4double cos_chi_exit;
  G4double phi_primary_cm;
  G4double cos_chi_secondary_8be;
  G4double cos_theta_secondary_correlation;
  G4double projectile_kinetic_lab;
  G4double projectile_px_lab;
  G4double projectile_py_lab;
  G4double projectile_pz_lab;
  G4double projectile_p_lab;
  G4double projectile_theta_lab;
  G4double projectile_phi_lab;

  G4int primary_angular_mode; // 0: isotropic, 1: fixed A1/A2
  G4int enable_162_primary_angular_distribution;
  G4double a1_162_primary;
  G4double a2_162_primary;
  G4int enable_675_primary_angular_distribution;
  G4double a1_675_primary;
  G4double a2_675_primary;
  G4double primary_a1_used;
  G4double primary_a2_used;

  G4int h11b675_decay_model; // kept for backward compatibility; 0: isotropic, 1: Legendre A2/A4, 2: strict coherent L1/L3
  G4int h11b675_decay_model_used;
  G4int h11b675_alpha_decay_model;
  G4int h11b675_strict_coherent_l13;
  G4int h11b675_strict_permutation_symmetrized;
  G4double h11b675_strict_l1_fraction;
  G4double h11b675_strict_l13_phase;
  G4double h11b675_strict_8be_lambda_energy_keV;
  G4double h11b675_strict_8be_reduced_width_squared_keV;
  G4double h11b675_strict_weight;
  G4double h11b675_strict_weight_max;
  G4int h11b675_strict_sampling_attempts;
  G4int enable_675_alpha1_secondary_angular_correlation;
  G4int secondary_angular_model; // 0: isotropic, 1: Legendre A2/A4
  G4double secondary_a2_used;
  G4double secondary_a4_used;
  G4int background_sequential_model; // kept for backward compatibility; current background is phase-space only
  G4double cos_chi_675_internal;
  G4double phi_chi_675_internal;
  G4double opening_angle_alpha12_cm;
  G4double opening_angle_alpha13_cm;
  G4double opening_angle_alpha23_cm;
  G4double e_alpha1_cm;
  G4double e_alpha2_cm;
  G4double e_alpha3_cm;
  G4double e_8be_excitation;

  G4double event_weight;

  // three-alpha center-of-mass kinetic energies
  // These are the quantities preferred for Dalitz plots.
  G4double e_3alpha_cm_alpha1;
  G4double e_3alpha_cm_alpha2;
  G4double e_3alpha_cm_alpha3;

  G4double ex_8Be; // exciation energy of 8Be
  G4double e_8Be;
  G4double theta_lab_8Be;
  G4double phi_lab_8Be;

  G4double x;
  G4double y;
  G4double z;

  G4double sigma_eval_b;
  G4double sigma_162_model_b;
  G4double sigma_675_model_b;
  G4double sigma_162_total_b;
  G4double sigma_675_total_b;
  G4double sigma_162_used_b;
  G4double sigma_675_used_b;
  G4double sigma_total_used_b;
  G4double sigma_162_sampling_b;
  G4double sigma_675_sampling_b;
  G4double sigma_162_directdecay_sampling_b;
  G4double sigma_675_directdecay_sampling_b;
  G4double sigma_background_sampling_b;
  G4double sigma_directdecay_sampling_b;
  G4double sigma_3alpha_sampling_total_b;
  G4double sigma_model_sum_b;
  G4double model_scale_factor;
  G4double cross_section_bias_factor;
  G4double background_bias_factor;
  G4double direct_decay_fraction;
  G4double sequential_decay_fraction_162;
  G4double sequential_decay_fraction_675;
  G4double direct_decay_fraction_162;
  G4double direct_decay_fraction_675;
  G4int enable_direct_decay;
  G4double scale_factor_162;
  G4double scale_factor_675;
  G4double sigma_background_b;
  G4double sigma_162_directdecay_b;
  G4double sigma_675_directdecay_b;
  G4double sigma_directdecay_b;
  G4double sigma_3alpha_eval_b;
  G4double sigma_gamma_162_0_b;
  G4double sigma_gamma_162_1_b;
  G4double sigma_gamma_162_total_b;
  G4double sigma_gamma_675_total_b;
  G4double sigma_gamma_total_b;
  G4double sigma_total_physical_all_b;
  G4double sigma_total_sampling_all_b;
  G4double sigma_gamma_162_0_physical_b;
  G4double sigma_gamma_162_1_physical_b;
  G4double sigma_gamma_675_physical_b;
  G4double sigma_gamma_162_0_sampling_b;
  G4double sigma_gamma_162_1_sampling_b;
  G4double sigma_gamma_675_sampling_b;
  G4double sigma_gamma_675_to_ground_physical_b;
  G4double sigma_gamma_675_to_4439_physical_b;
  G4double sigma_gamma_675_to_7654_physical_b;
  G4double sigma_gamma_675_to_12710_physical_b;
  G4double sigma_gamma_675_to_15110_physical_b;
  G4double sigma_gamma_675_to_ground_sampling_b;
  G4double sigma_gamma_675_to_4439_sampling_b;
  G4double sigma_gamma_675_to_7654_sampling_b;
  G4double sigma_gamma_675_to_12710_sampling_b;
  G4double sigma_gamma_675_to_15110_sampling_b;
  G4double sigma_gamma_sampling_total_b;
  G4double sigma_total_all_b;
  G4double channel_probability_162;
  G4double channel_probability_675;
  G4double channel_probability_background;
  G4double channel_probability_directdecay;
  G4double channel_probability_directdecay_162;
  G4double channel_probability_directdecay_675;
  G4double channel_probability_gamma;
  G4double probability_gamma_162_0;
  G4double probability_gamma_162_1;
  G4double probability_gamma_675_total;
  G4double probability_gamma_675_to_ground;
  G4double probability_gamma_675_to_4439;
  G4double probability_gamma_675_to_7654;
  G4double probability_gamma_675_to_12710;
  G4double probability_gamma_675_to_15110;

  // phaseSpace
  char reaction[128];

  //
  H11BReactionData()
      : event(0), e_alpha1(0), e_alpha2(0), e_alpha3(0), theta_lab_alpha1(0), theta_lab_alpha2(0), theta_lab_alpha3(0),
        phi_lab_alpha1(0), phi_lab_alpha2(0), phi_lab_alpha3(0), resonance_id(0), branch_id(0),
        reaction_channel(0), background_mode(0), gamma_resonance(0), gamma_branch(0), gamma_angular_mode(0),
        n_prompt_gammas(0), gamma1_energy(0), gamma2_energy(0), gamma1_theta_lab(0), gamma2_theta_lab(0),
        gamma1_phi_lab(0), gamma2_phi_lab(0), gamma1_theta_cm(0), gamma2_theta_cm(0), cos_theta_gamma_cm(0),
        gamma_event_weight(1), gamma_bias_factor(1), event_sampling_weight(1), gamma_final_state_energy_MeV(0),
        gamma_primary_energy_MeV(0), gamma_relative_intensity_used(0), gamma_branch_fraction_used(0),
        gamma_branch_is_upper_limit(0), gamma_branch_from_relative_table(0), gamma_cascade_generated(0),
        enable_675_gamma_angular_distribution(0), a1_675_gamma(0), a2_675_gamma(0),
        e_cm_p11B(0), ex_max_8Be(0), eaa_8Be(0), e_alpha8Be(0),
        cos_theta_primary_cm(0), cos_chi_exit(0), phi_primary_cm(0), cos_chi_secondary_8be(0),
        cos_theta_secondary_correlation(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_kinetic_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_px_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_py_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_pz_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_p_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_theta_lab(std::numeric_limits<G4double>::quiet_NaN()),
        projectile_phi_lab(std::numeric_limits<G4double>::quiet_NaN()),
        primary_angular_mode(0), enable_162_primary_angular_distribution(0), a1_162_primary(0), a2_162_primary(0),
        enable_675_primary_angular_distribution(0), a1_675_primary(0), a2_675_primary(0),
        primary_a1_used(0), primary_a2_used(0), h11b675_decay_model(0),
        h11b675_decay_model_used(0), h11b675_alpha_decay_model(0), h11b675_strict_coherent_l13(0),
        h11b675_strict_permutation_symmetrized(0), h11b675_strict_l1_fraction(0), h11b675_strict_l13_phase(0),
        h11b675_strict_8be_lambda_energy_keV(0), h11b675_strict_8be_reduced_width_squared_keV(0),
        h11b675_strict_weight(0), h11b675_strict_weight_max(0), h11b675_strict_sampling_attempts(0),
        enable_675_alpha1_secondary_angular_correlation(0), secondary_angular_model(0), secondary_a2_used(0), secondary_a4_used(0),
        background_sequential_model(0), cos_chi_675_internal(0),
        phi_chi_675_internal(0), opening_angle_alpha12_cm(0), opening_angle_alpha13_cm(0),
        opening_angle_alpha23_cm(0), e_alpha1_cm(0), e_alpha2_cm(0), e_alpha3_cm(0), e_8be_excitation(0),
        event_weight(1), e_3alpha_cm_alpha1(0), e_3alpha_cm_alpha2(0), e_3alpha_cm_alpha3(0), ex_8Be(0),
        e_8Be(0), theta_lab_8Be(0), phi_lab_8Be(0),

        x(0.), y(0.), z(0.), sigma_eval_b(0.), sigma_162_model_b(0.), sigma_675_model_b(0.),
        sigma_162_total_b(0.), sigma_675_total_b(0.),
        sigma_162_used_b(0.), sigma_675_used_b(0.), sigma_total_used_b(0.),
        sigma_162_sampling_b(0.), sigma_675_sampling_b(0.), sigma_162_directdecay_sampling_b(0.),
        sigma_675_directdecay_sampling_b(0.), sigma_background_sampling_b(0.),
        sigma_directdecay_sampling_b(0.), sigma_3alpha_sampling_total_b(0.), sigma_model_sum_b(0.), model_scale_factor(1.),
        cross_section_bias_factor(1.), background_bias_factor(0.01), direct_decay_fraction(0.01),
        sequential_decay_fraction_162(0.99), sequential_decay_fraction_675(0.99),
        direct_decay_fraction_162(0.01), direct_decay_fraction_675(0.01), enable_direct_decay(1),
        scale_factor_162(1.), scale_factor_675(1.), sigma_background_b(0.), sigma_162_directdecay_b(0.),
        sigma_675_directdecay_b(0.), sigma_directdecay_b(0.),
        sigma_3alpha_eval_b(0.), sigma_gamma_162_0_b(0.),
        sigma_gamma_162_1_b(0.), sigma_gamma_162_total_b(0.), sigma_gamma_675_total_b(0.),
        sigma_gamma_total_b(0.), sigma_total_physical_all_b(0.), sigma_total_sampling_all_b(0.),
        sigma_gamma_162_0_physical_b(0.), sigma_gamma_162_1_physical_b(0.), sigma_gamma_675_physical_b(0.),
        sigma_gamma_162_0_sampling_b(0.), sigma_gamma_162_1_sampling_b(0.), sigma_gamma_675_sampling_b(0.),
        sigma_gamma_675_to_ground_physical_b(0.), sigma_gamma_675_to_4439_physical_b(0.),
        sigma_gamma_675_to_7654_physical_b(0.), sigma_gamma_675_to_12710_physical_b(0.),
        sigma_gamma_675_to_15110_physical_b(0.), sigma_gamma_675_to_ground_sampling_b(0.),
        sigma_gamma_675_to_4439_sampling_b(0.), sigma_gamma_675_to_7654_sampling_b(0.),
        sigma_gamma_675_to_12710_sampling_b(0.), sigma_gamma_675_to_15110_sampling_b(0.),
        sigma_gamma_sampling_total_b(0.), sigma_total_all_b(0.),
        channel_probability_162(0.), channel_probability_675(0.), channel_probability_background(0.),
        channel_probability_directdecay(0.), channel_probability_directdecay_162(0.),
        channel_probability_directdecay_675(0.), channel_probability_gamma(0.), probability_gamma_162_0(0.), probability_gamma_162_1(0.),
        probability_gamma_675_total(0.), probability_gamma_675_to_ground(0.), probability_gamma_675_to_4439(0.),
        probability_gamma_675_to_7654(0.), probability_gamma_675_to_12710(0.), probability_gamma_675_to_15110(0.)
  {
    std::memset(reaction, '\0', sizeof(reaction));
  }

  //
  void Clear()
  {
    *this = H11BReactionData();
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct EventData
{
  G4long event_id;

  // Independent DSSD strip signals.  One vector element is one electronic
  // channel after all energy deposits on that strip have been accumulated.
  std::vector<G4int> si_detector_id;
  std::vector<G4int> si_side;     // 0: front, 1: back
  std::vector<G4int> si_strip_id;
  std::vector<G4double> si_edep_MeV;
  std::vector<G4double> si_time_ns;

  // Unsegmented gamma-detector channel signals.
  std::vector<G4int> labr3_detector_id;
  std::vector<G4double> labr3_edep_MeV;
  std::vector<G4double> labr3_time_ns;

  std::vector<G4int> hpge_detector_id;
  std::vector<G4double> hpge_edep_MeV;
  std::vector<G4double> hpge_time_ns;

  EventData()
      : event_id(-1)
  {
  }

  void Clear()
  {
    event_id = -1;

    si_detector_id.clear();
    si_side.clear();
    si_strip_id.clear();
    si_edep_MeV.clear();
    si_time_ns.clear();

    labr3_detector_id.clear();
    labr3_edep_MeV.clear();
    labr3_time_ns.clear();

    hpge_detector_id.clear();
    hpge_edep_MeV.clear();
    hpge_time_ns.clear();
  }
};

// Static DSSD pixel-centre geometry.  These entries contain detector geometry
// only; they do not contain event or Monte-Carlo track truth.
struct SiPixelMapData
{
  G4int detector_id;
  G4int detector_model; // 1: W1, 2: S3
  G4int subarray_id;
  G4int module_id;
  G4int front_strip_id;
  G4int back_strip_id;

  G4double x_center_mm;
  G4double y_center_mm;
  G4double z_center_mm;
  G4double theta_lab_center_deg;
  G4double phi_lab_center_deg;

  SiPixelMapData()
      : detector_id(-1), detector_model(0), subarray_id(-1), module_id(-1), front_strip_id(-1),
        back_strip_id(-1), x_center_mm(0.), y_center_mm(0.), z_center_mm(0.), theta_lab_center_deg(0.),
        phi_lab_center_deg(0.)
  {
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct TrackData
{
  G4long event;
  G4int track;
  G4double e;
  G4double x;
  G4double y;
  G4double z;
  G4double ts;
  G4double length;
  char volume[128];
  char particle[128];

  //
  TrackData()
      : event(-1), track(-1), e(0.), x(0.), y(0.), z(0.), ts(0.), length(0.)
  {
    std::memset(volume, '\0', sizeof(volume));
    std::memset(particle, '\0', sizeof(particle));
  }

  //
  void Clear()
  {
    *this = TrackData();
  }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct StepData
{
  G4long event;
  G4int track;
  G4double de;
  G4double pre_x;
  G4double pre_y;
  G4double pre_z;
  G4double pre_total_energy;
  G4double pre_kine_energy;
  G4double post_x;
  G4double post_y;
  G4double post_z;
  G4double post_total_energy;
  G4double post_kine_energy;
  G4double length;
  char volume[128];
  char particle[128];
  char process[128];

  //
  StepData()
      : event(-1), track(-1), de(0.), pre_x(0.), pre_y(0.), pre_z(0.), pre_total_energy(0.), pre_kine_energy(0.),
        post_x(0.), post_y(0.), post_z(0.), post_total_energy(0.), post_kine_energy(0.), length(0.)
  {
    std::memset(volume, '\0', sizeof(volume));
    std::memset(particle, '\0', sizeof(particle));
    std::memset(process, '\0', sizeof(process));
  }

  //
  void Clear()
  {
    *this = StepData();
  }
};
#endif
