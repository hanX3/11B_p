#include "H11BCrossSection.hh"

#include "Constants.hh"
#include "CoulombPenetrability.hh"
#include "H11BConfig.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {
G4double g_cross_section_bias_factor = 1.0;
G4double g_162_sequential_decay_fraction = 0.99;
G4double g_675_sequential_decay_fraction = 0.99;
G4bool g_enable_direct_decay = true;
G4double g_162_bw_scale_factor = 1.0;
G4double g_675_scale_factor = 1.0;

constexpr G4double H11B675GammaRelativeIntensitySum = 15.7 + 100.0 + 6.8 + 0.16;

G4double ClampBranchFraction(G4double fraction)
{
  return std::clamp(fraction, 0.0, 1.0);
}

void Fill675GammaPhysicalBranches(H11BCrossSectionComponents& components)
{
  if (components.sigma_gamma_675_physical <= 0. || H11B675GammaRelativeIntensitySum <= 0.)
    return;

  components.sigma_gamma_675_to_ground_physical = components.sigma_gamma_675_physical * 15.7 / H11B675GammaRelativeIntensitySum;
  components.sigma_gamma_675_to_4439_physical = components.sigma_gamma_675_physical * 100.0 / H11B675GammaRelativeIntensitySum;
  components.sigma_gamma_675_to_12710_physical = components.sigma_gamma_675_physical * 6.8 / H11B675GammaRelativeIntensitySum;
  components.sigma_gamma_675_to_15110_physical = components.sigma_gamma_675_physical * 0.16 / H11B675GammaRelativeIntensitySum;

  if (H11B675IncludeUpperLimitGammaLines) {
    // Keep the upper-limit branch out of the default normalization.  If enabled
    // for systematic checks, assign its small intensity relative to the same
    // default sum so the primary branches remain unchanged.
    components.sigma_gamma_675_to_7654_physical = components.sigma_gamma_675_physical * 0.07 / H11B675GammaRelativeIntensitySum;
  }
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool H11BCrossSection::IsIsoApplicable(const G4DynamicParticle*, G4int, G4int, const G4Element*, const G4Material*)
{
  // This is the only data set attached to the custom proton process, so it
  // must answer every isotope query. Non-11B targets receive zero below.
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetIsoCrossSection(const G4DynamicParticle* projectile, G4int z, G4int a, const G4Isotope*, const G4Element*, const G4Material*)
{
  if (!projectile || projectile->GetDefinition() != G4Proton::Definition() || z != 5 || a != 11) return 0.;
  return CalculateComponents(projectile->GetKineticEnergy()).sigma_total_sampling_all;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H11BCrossSectionComponents H11BCrossSection::CalculateComponents(G4double kinetic_energy_lab)
{
  const G4double energy_cm_keV = GetEcmValue(1., 11., kinetic_energy_lab / keV);

  H11BCrossSectionComponents components;
  components.sigma_162_model = GetSigma162(energy_cm_keV) * cm2;

  // Phenomenological weighted_chebE16 fit675 component from cs_model.  This is
  // not an isolated 675-keV Breit-Wigner resonance.
  components.sigma_675_model = GetFit675CrossSection(kinetic_energy_lab);

  components.sigma_162_total = std::max(0.0, g_162_bw_scale_factor) * components.sigma_162_model;
  components.sigma_675_total = std::max(0.0, g_675_scale_factor) * components.sigma_675_model;

  const G4double seq162 = g_enable_direct_decay ? ClampBranchFraction(g_162_sequential_decay_fraction) : 1.0;
  const G4double seq675 = g_enable_direct_decay ? ClampBranchFraction(g_675_sequential_decay_fraction) : 1.0;
  const G4double direct162 = g_enable_direct_decay ? 1.0 - seq162 : 0.0;
  const G4double direct675 = g_enable_direct_decay ? 1.0 - seq675 : 0.0;

  components.sigma_162 = seq162 * components.sigma_162_total;
  components.sigma_675 = seq675 * components.sigma_675_total;
  components.sigma_162_directdecay = direct162 * components.sigma_162_total;
  components.sigma_675_directdecay = direct675 * components.sigma_675_total;
  components.sigma_directdecay = components.sigma_162_directdecay + components.sigma_675_directdecay;

  // Backward-compatible aliases for old analysis scripts.
  components.sigma_background = components.sigma_directdecay;

  components.sigma_total = components.sigma_162 + components.sigma_675 + components.sigma_directdecay;
  // Backward-compatible diagnostic alias for the only remaining 3-alpha model:
  // scaled 162 BW plus scaled fit675, independent of the sequential/direct split.
  components.sigma_eval = components.sigma_total;
  components.cross_section_bias_factor = g_cross_section_bias_factor;
  components.direct_decay_fraction =
    components.sigma_total > 0.0 ? components.sigma_directdecay / components.sigma_total : 0.0;
  components.sequential_decay_fraction_162 = seq162;
  components.sequential_decay_fraction_675 = seq675;
  components.direct_decay_fraction_162 = direct162;
  components.direct_decay_fraction_675 = direct675;
  components.scale_factor_162 = g_162_bw_scale_factor;
  components.scale_factor_675 = g_675_scale_factor;
  components.enable_direct_decay = g_enable_direct_decay;
  components.background_bias_factor = components.direct_decay_fraction;

  components.sigma_162_sampling = components.sigma_162 * components.cross_section_bias_factor;
  components.sigma_675_sampling = components.sigma_675 * components.cross_section_bias_factor;
  components.sigma_162_directdecay_sampling = components.sigma_162_directdecay * components.cross_section_bias_factor;
  components.sigma_675_directdecay_sampling = components.sigma_675_directdecay * components.cross_section_bias_factor;
  components.sigma_directdecay_sampling = components.sigma_directdecay * components.cross_section_bias_factor;
  components.sigma_background_sampling = components.sigma_directdecay_sampling;
  components.sigma_3alpha_sampling_total =
    components.sigma_162_sampling + components.sigma_675_sampling + components.sigma_directdecay_sampling;

  if (H11BConfig::Get162GammaCaptureEnabled()) {
    components.sigma_gamma_162_0_physical = Get162Gamma0CrossSection(energy_cm_keV) * cm2;
    components.sigma_gamma_162_1_physical = Get162Gamma1CrossSection(energy_cm_keV) * cm2;
    components.sigma_gamma_162_total = components.sigma_gamma_162_0_physical + components.sigma_gamma_162_1_physical;
  }

  if (H11BConfig::Get675GammaCaptureEnabled()) {
    // Gamma capture is an independent exit channel (README "Gamma Capture
    // Channels": sigma_gamma_675 = 1e-5 * sigma_675_model), so it must use
    // sigma_675_model, not the sequential-decay-fraction-scaled sigma_675 --
    // otherwise 675SequentialDecayFraction=0 zeroes out gamma capture too.
    components.sigma_gamma_675_physical = Get675GammaTotalCrossSection(components.sigma_675_model);
    Fill675GammaPhysicalBranches(components);
  }

  components.sigma_gamma_total_physical = components.sigma_gamma_162_total + components.sigma_gamma_675_physical;

  // Backward-compatible aliases.  These are physical cross sections.
  components.sigma_gamma_162_0 = components.sigma_gamma_162_0_physical;
  components.sigma_gamma_162_1 = components.sigma_gamma_162_1_physical;
  components.sigma_gamma_675_total = components.sigma_gamma_675_physical;
  components.sigma_gamma_total = components.sigma_gamma_total_physical;

  const G4double bias_factor = H11BConfig::GetGammaBiasFactor();
  const G4double total_bias_factor = bias_factor * components.cross_section_bias_factor;
  components.sigma_gamma_162_0_sampling = components.sigma_gamma_162_0_physical * total_bias_factor;
  components.sigma_gamma_162_1_sampling = components.sigma_gamma_162_1_physical * total_bias_factor;
  components.sigma_gamma_675_sampling = components.sigma_gamma_675_physical * total_bias_factor;
  components.sigma_gamma_675_to_ground_sampling = components.sigma_gamma_675_to_ground_physical * total_bias_factor;
  components.sigma_gamma_675_to_4439_sampling = components.sigma_gamma_675_to_4439_physical * total_bias_factor;
  components.sigma_gamma_675_to_7654_sampling = components.sigma_gamma_675_to_7654_physical * total_bias_factor;
  components.sigma_gamma_675_to_12710_sampling = components.sigma_gamma_675_to_12710_physical * total_bias_factor;
  components.sigma_gamma_675_to_15110_sampling = components.sigma_gamma_675_to_15110_physical * total_bias_factor;
  components.sigma_gamma_sampling_total = components.sigma_gamma_162_0_sampling + components.sigma_gamma_162_1_sampling + components.sigma_gamma_675_sampling;

  components.sigma_total_physical_all = components.sigma_total + components.sigma_gamma_total_physical;
  components.sigma_total_sampling_all = components.sigma_3alpha_sampling_total + components.sigma_gamma_sampling_total;
  components.sigma_total_all = components.sigma_total_sampling_all;

  return components;
}

G4double H11BCrossSection::GetCrossSectionBiasFactor()
{
  return g_cross_section_bias_factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::SetCrossSectionBiasFactor(G4double factor)
{
  if (!std::isfinite(factor) || factor < 1.0) {
    G4cerr << "Invalid /h11b/crossSectionBiasFactor " << factor << ". Use a finite value >= 1." << G4endl;
    return;
  }

  g_cross_section_bias_factor = factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get162SequentialDecayFraction()
{
  return g_162_sequential_decay_fraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::Set162SequentialDecayFraction(G4double fraction)
{
  if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
    G4cerr << "Invalid /h11b/162SequentialDecayFraction " << fraction << ". Use a finite value in [0, 1]." << G4endl;
    return;
  }

  g_162_sequential_decay_fraction = fraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get675SequentialDecayFraction()
{
  return g_675_sequential_decay_fraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::Set675SequentialDecayFraction(G4double fraction)
{
  if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
    G4cerr << "Invalid /h11b/675SequentialDecayFraction " << fraction << ". Use a finite value in [0, 1]." << G4endl;
    return;
  }

  g_675_sequential_decay_fraction = fraction;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool H11BCrossSection::GetDirectDecayEnabled()
{
  return g_enable_direct_decay;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::SetDirectDecayEnabled(G4bool enabled)
{
  g_enable_direct_decay = enabled;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get162BWScaleFactor()
{
  return g_162_bw_scale_factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::Set162BWScaleFactor(G4double factor)
{
  if (!std::isfinite(factor) || factor < 0.0) {
    G4cerr << "Invalid /h11b/162BWScaleFactor " << factor << ". Use a finite value >= 0." << G4endl;
    return;
  }

  g_162_bw_scale_factor = factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get675ScaleFactor()
{
  return g_675_scale_factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::Set675ScaleFactor(G4double factor)
{
  if (!std::isfinite(factor) || factor < 0.0) {
    G4cerr << "Invalid /h11b/675ScaleFactor " << factor << ". Use a finite value >= 0." << G4endl;
    return;
  }

  g_675_scale_factor = factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetEcmValue(G4double project_a, G4double target_a, G4double kinetic_energy_lab_keV)
{
  return target_a / (project_a + target_a) * kinetic_energy_lab_keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetSigma162(G4double ecm)
{
  if (ecm < 1.) return 0.;

  G4double exit_width = H11B162Alpha0Width / keV + H11B162Alpha1Width / keV;
  if (H11BIncludeGammaChannelInCrossSection) {
    exit_width += H11B162GammaWidth / keV;
  }

  return GetSigmaBreitWigner(ecm, H11B162ResonanceEnergy / keV, H11B162TotalWidth / keV, H11B162ProtonWidth / keV, exit_width, H11B162SpinStatFactor, H11B162EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get162Gamma0CrossSection(G4double ecm)
{
  if (ecm < 1. || ecm > 400.) return 0.;

  return GetSigmaBreitWigner(ecm, H11B162ResonanceEnergy / keV, H11B162TotalWidth / keV, H11B162ProtonWidth / keV, H11B162Gamma0Width / keV, H11B162SpinStatFactor, H11B162EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get162Gamma1CrossSection(G4double ecm)
{
  if (ecm < 1. || ecm > 400.) return 0.;

  return GetSigmaBreitWigner(ecm, H11B162ResonanceEnergy / keV, H11B162TotalWidth / keV, H11B162ProtonWidth / keV, H11B162Gamma1Width / keV, H11B162SpinStatFactor, H11B162EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get675GammaTotalCrossSection(G4double sigma_675_effective)
{
  if (!H11BConfig::Get675GammaCaptureEnabled() || sigma_675_effective <= 0.) return 0.;

  return H11B675GammaToAlphaScale * sigma_675_effective;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetFit675CrossSection(G4double kinetic_energy_lab)
{
  if (kinetic_energy_lab <= 0.) return 0.;

  constexpr G4double fit_min_lab_MeV = 0.020;
  constexpr G4double fit_max_lab_MeV = 1.000;
  constexpr G4double coefficients[] = {
    -4.240991093218323,
    6.5335030639422,
    -5.348623593275032,
    2.836731408213122,
    -1.88972601027126,
    1.556619927270825,
    -0.9501800114612392,
    0.6222767937946161,
    -0.4865823078410484,
    0.3336851443104574,
    -0.2357424108239984,
    0.1936547242371347,
    -0.1385260429805155,
    0.1008831241333581,
    -0.05830063682899136,
    0.02350578311364693,
    0.009138850080015581,
  };

  // Analytic weighted_chebE16 fit675 from cs_model.  The fit is calibrated on
  // 0.020-1.000 MeV proton lab energy; clamp outside that range to avoid
  // uncontrolled Chebyshev extrapolation in higher-energy runs.
  const G4double e_lab_MeV = std::clamp(kinetic_energy_lab / MeV, fit_min_lab_MeV, fit_max_lab_MeV);
  const G4double z = 2.0 * (e_lab_MeV - fit_min_lab_MeV) / (fit_max_lab_MeV - fit_min_lab_MeV) - 1.0;

  G4double log_sigma_b = coefficients[0];
  G4double t_nm2 = 1.0;
  G4double t_nm1 = z;
  log_sigma_b += coefficients[1] * t_nm1;
  for (std::size_t i = 2; i < sizeof(coefficients) / sizeof(coefficients[0]); ++i) {
    const G4double t_n = 2.0 * z * t_nm1 - t_nm2;
    log_sigma_b += coefficients[i] * t_n;
    t_nm2 = t_nm1;
    t_nm1 = t_n;
  }

  return std::exp(log_sigma_b) * barn;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetSigmaBreitWigner(G4double ecm, G4double resonance_energy, G4double total_width, G4double entrance_width_at_resonance, G4double exit_width_at_resonance, G4double spin_stat_factor, G4int entrance_orbital_l)
{
  if (ecm <= 0. || resonance_energy <= 0. || total_width <= 0. || entrance_width_at_resonance <= 0. || exit_width_at_resonance <= 0.) return 0.;

  // Single-level Breit-Wigner formula with an energy-dependent p+11B
  // entrance width:
  //
  //   Gamma_p(E) = Gamma_p(Er) * P_L(E,a) / P_L(Er,a)
  //
  // where P_L is the full Coulomb penetrability
  //
  //   P_L = rho / [F_L(eta,rho)^2 + G_L(eta,rho)^2].
  //
  // The ratio is evaluated from precomputed Coulomb-function tables.  The
  // remaining non-entrance width is kept constant, so the total width at Er
  // is unchanged.
  const G4double projectA = 1.;
  const G4double targetA = 11.;
  const G4double amu_c2_keV = 931494.10242;
  const G4double hbarc_keV_fm = 197326.9804;
  const G4double reduced_mass_keV = projectA * targetA / (projectA + targetA) * amu_c2_keV;

  const G4double pi_over_k2_fm2 = pi * hbarc_keV_fm * hbarc_keV_fm / (2. * reduced_mass_keV * ecm);

  const G4double entrance_width_e = GetFullCoulombEntranceWidth(ecm, resonance_energy, entrance_width_at_resonance, entrance_orbital_l);

  const G4double missing_width = std::max(0., total_width - entrance_width_at_resonance - exit_width_at_resonance);
  const G4double total_width_e = entrance_width_e + exit_width_at_resonance + missing_width;

  const G4double denominator = (ecm - resonance_energy) * (ecm - resonance_energy) + total_width_e * total_width_e / 4.;
  const G4double bw_factor = spin_stat_factor * entrance_width_e * exit_width_at_resonance / denominator;

  return pi_over_k2_fm2 * bw_factor * 1e-26; // fm2 -> cm2
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetFullCoulombEntranceWidth(G4double ecm, G4double resonance_energy, G4double entrance_width_at_resonance, G4int entrance_orbital_l)
{
  if (ecm <= 0. || resonance_energy <= 0. || entrance_width_at_resonance <= 0.) return 0.;

  const G4double ratio = CoulombPenetrability::P11BRatio(ecm, resonance_energy, entrance_orbital_l);

  return entrance_width_at_resonance * ratio;
}
