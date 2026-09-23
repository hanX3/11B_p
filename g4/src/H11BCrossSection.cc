#include "H11BCrossSection.hh"

#include "Constants.hh"
#include "CoulombPenetrability.hh"
#include "H11BConfig.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cmath>

namespace {
H11BEvaluatedCrossSectionMode g_evaluated_cross_section_mode = H11BEvaluatedCrossSectionMode::Tentori2023;
G4double g_cross_section_bias_factor = 1.0;
G4double g_background_bias_factor = 1.0;

constexpr G4double LabToCm = 11.0 / 12.0;
constexpr G4double GamowEnergyMeV = 22.589;
constexpr G4double H11B675GammaRelativeIntensitySum = 15.7 + 100.0 + 6.8 + 0.16;

G4double Square(G4double x)
{
  return x * x;
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
  components.sigma_165_model = GetSigma165(energy_cm_keV) * cm2;
  components.sigma_675_model = GetSigma675(energy_cm_keV) * cm2;

  const G4double sigma_model_sum = components.sigma_165_model + components.sigma_675_model;
  components.sigma_165 = components.sigma_165_model;
  components.sigma_675 = components.sigma_675_model;

  if (g_evaluated_cross_section_mode == H11BEvaluatedCrossSectionMode::Model) {
    components.sigma_eval = sigma_model_sum;
  } else {
    components.sigma_eval = GetEvaluatedTotalCrossSection(kinetic_energy_lab);
  }

  if (components.sigma_eval <= 0.) {
    components.sigma_165 = 0.;
    components.sigma_675 = 0.;
    components.sigma_background = 0.;
  } else if (sigma_model_sum <= components.sigma_eval) {
    components.sigma_background = components.sigma_eval - sigma_model_sum;
  } else {
    const G4double scale = components.sigma_eval / sigma_model_sum;
    components.sigma_165 *= scale;
    components.sigma_675 *= scale;
    components.sigma_background = 0.;
  }

  components.sigma_total = components.sigma_165 + components.sigma_675 + components.sigma_background;
  components.cross_section_bias_factor = g_cross_section_bias_factor;
  components.background_bias_factor = g_background_bias_factor;
  components.sigma_165_sampling = components.sigma_165 * components.cross_section_bias_factor;
  components.sigma_675_sampling = components.sigma_675 * components.cross_section_bias_factor;
  components.sigma_background_sampling =
    components.sigma_background * components.cross_section_bias_factor * components.background_bias_factor;
  components.sigma_3alpha_sampling_total =
    components.sigma_165_sampling + components.sigma_675_sampling + components.sigma_background_sampling;

  if (H11BConfig::Get165GammaCaptureEnabled()) {
    components.sigma_gamma_165_0_physical = Get165Gamma0CrossSection(energy_cm_keV) * cm2;
    components.sigma_gamma_165_1_physical = Get165Gamma1CrossSection(energy_cm_keV) * cm2;
    components.sigma_gamma_165_total = components.sigma_gamma_165_0_physical + components.sigma_gamma_165_1_physical;
  }

  if (H11BConfig::Get675GammaCaptureEnabled()) {
    components.sigma_gamma_675_physical = Get675GammaTotalCrossSection(components.sigma_675_model);
    Fill675GammaPhysicalBranches(components);
  }

  components.sigma_gamma_total_physical = components.sigma_gamma_165_total + components.sigma_gamma_675_physical;

  // Backward-compatible aliases.  These are physical cross sections.
  components.sigma_gamma_165_0 = components.sigma_gamma_165_0_physical;
  components.sigma_gamma_165_1 = components.sigma_gamma_165_1_physical;
  components.sigma_gamma_675_total = components.sigma_gamma_675_physical;
  components.sigma_gamma_total = components.sigma_gamma_total_physical;

  const G4double bias_factor = H11BConfig::GetGammaBiasFactor();
  const G4double total_bias_factor = bias_factor * components.cross_section_bias_factor;
  components.sigma_gamma_165_0_sampling = components.sigma_gamma_165_0_physical * total_bias_factor;
  components.sigma_gamma_165_1_sampling = components.sigma_gamma_165_1_physical * total_bias_factor;
  components.sigma_gamma_675_sampling = components.sigma_gamma_675_physical * total_bias_factor;
  components.sigma_gamma_675_to_ground_sampling = components.sigma_gamma_675_to_ground_physical * total_bias_factor;
  components.sigma_gamma_675_to_4439_sampling = components.sigma_gamma_675_to_4439_physical * total_bias_factor;
  components.sigma_gamma_675_to_7654_sampling = components.sigma_gamma_675_to_7654_physical * total_bias_factor;
  components.sigma_gamma_675_to_12710_sampling = components.sigma_gamma_675_to_12710_physical * total_bias_factor;
  components.sigma_gamma_675_to_15110_sampling = components.sigma_gamma_675_to_15110_physical * total_bias_factor;
  components.sigma_gamma_sampling_total = components.sigma_gamma_165_0_sampling + components.sigma_gamma_165_1_sampling + components.sigma_gamma_675_sampling;

  components.sigma_total_physical_all = components.sigma_total + components.sigma_gamma_total_physical;
  components.sigma_total_sampling_all = components.sigma_3alpha_sampling_total + components.sigma_gamma_sampling_total;
  components.sigma_total_all = components.sigma_total_sampling_all;

  return components;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H11BEvaluatedCrossSectionMode H11BCrossSection::GetEvaluatedCrossSectionMode()
{
  return g_evaluated_cross_section_mode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::SetEvaluatedCrossSectionMode(H11BEvaluatedCrossSectionMode mode)
{
  g_evaluated_cross_section_mode = mode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::SetEvaluatedCrossSectionMode(const G4String& mode)
{
  if (mode == "tentori2023") {
    SetEvaluatedCrossSectionMode(H11BEvaluatedCrossSectionMode::Tentori2023);
  } else if (mode == "model") {
    SetEvaluatedCrossSectionMode(H11BEvaluatedCrossSectionMode::Model);
  } else {
    G4cerr << "Unknown /h11b/evaluatedCrossSection mode '" << mode << "'. Use tentori2023 or model." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const char* H11BCrossSection::GetEvaluatedCrossSectionModeName()
{
  switch (g_evaluated_cross_section_mode) {
  case H11BEvaluatedCrossSectionMode::Tentori2023:
    return "tentori2023";
  case H11BEvaluatedCrossSectionMode::Model:
    return "model";
  }

  return "unknown";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
G4double H11BCrossSection::GetBackgroundBiasFactor()
{
  return g_background_bias_factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BCrossSection::SetBackgroundBiasFactor(G4double factor)
{
  if (!std::isfinite(factor) || factor < 0.0) {
    G4cerr << "Invalid /h11b/backgroundBiasFactor " << factor << ". Use a finite value >= 0." << G4endl;
    return;
  }

  g_background_bias_factor = factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetEcmValue(G4double project_a, G4double target_a, G4double kinetic_energy_lab_keV)
{
  return target_a / (project_a + target_a) * kinetic_energy_lab_keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetSigma165(G4double ecm)
{
  if (ecm < 1. || ecm > 400.) return 0.;

  G4double exit_width = H11B165Alpha0Width / keV + H11B165Alpha1Width / keV;
  if (H11BIncludeGammaChannelInCrossSection) {
    exit_width += H11B165GammaWidth / keV;
  }

  return GetSigmaBreitWigner(ecm, H11B165ResonanceEnergy / keV, H11B165TotalWidth / keV, H11B165ProtonWidth / keV, exit_width, H11B165SpinStatFactor, H11B165EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetSigma675(G4double ecm)
{
  if (ecm < 1. || ecm > 3500.) return 0.;

  G4double exit_width = H11B675Alpha0Width / keV + H11B675Alpha1Width / keV;
  if (H11BIncludeGammaChannelInCrossSection) {
    exit_width += (H11B675Gamma0WidthReferenceUnused + H11B675Gamma1WidthReferenceUnused) / keV;
  }

  return GetSigmaBreitWigner(ecm, H11B675ResonanceEnergy / keV, H11B675TotalWidth / keV, H11B675ProtonWidth / keV, exit_width, H11B675SpinStatFactor, H11B675EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get165Gamma0CrossSection(G4double ecm)
{
  if (ecm < 1. || ecm > 400.) return 0.;

  return GetSigmaBreitWigner(ecm, H11B165ResonanceEnergy / keV, H11B165TotalWidth / keV, H11B165ProtonWidth / keV, H11B165Gamma0Width / keV, H11B165SpinStatFactor, H11B165EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get165Gamma1CrossSection(G4double ecm)
{
  if (ecm < 1. || ecm > 400.) return 0.;

  return GetSigmaBreitWigner(ecm, H11B165ResonanceEnergy / keV, H11B165TotalWidth / keV, H11B165ProtonWidth / keV, H11B165Gamma1Width / keV, H11B165SpinStatFactor, H11B165EntranceOrbitalL);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::Get675GammaTotalCrossSection(G4double sigma_675_model)
{
  if (!H11BConfig::Get675GammaCaptureEnabled() || sigma_675_model <= 0.) return 0.;

  return H11B675GammaToAlphaScale * sigma_675_model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetEvaluatedTotalCrossSection(G4double kinetic_energy_lab)
{
  const G4double energy_cm_MeV = LabToCm * kinetic_energy_lab / MeV;
  if (energy_cm_MeV <= 0.) return 0.;

  return GetTentori2023TotalCrossSection(energy_cm_MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetTentori2023TotalCrossSection(G4double energy_cm_MeV)
{
  const G4double e = energy_cm_MeV;
  const G4double x_keV = 1000.0 * e;
  G4double s_factor = 0.;

  if (e <= 0.400) {
    constexpr G4double C0 = 197.0;
    constexpr G4double C1 = 0.269;
    constexpr G4double C2 = 2.54e-4;
    constexpr G4double AL = 1.82e4;
    constexpr G4double EL = 148.0;
    constexpr G4double dEL = 2.35;

    s_factor = C0 + C1 * x_keV + C2 * x_keV * x_keV;
    s_factor += AL / (Square(x_keV - EL) + Square(dEL));
  } else if (e <= 0.668) {
    constexpr G4double D0 = 346.0;
    constexpr G4double D1 = 150.0;
    constexpr G4double D2 = -59.9;
    constexpr G4double D5 = -0.460;

    const G4double x = (x_keV - 400.0) / 100.0;
    s_factor = D0 + D1 * x + D2 * x * x + D5 * std::pow(x, 5.0);
  } else {
    constexpr G4double B = 0.381;
    constexpr G4double A[] = {1.98e6, 3.89e6, 1.36e6, 3.71e6};
    constexpr G4double ER[] = {640.9, 1211.0, 2340.0, 3294.0};
    constexpr G4double dE[] = {85.5, 414.0, 221.0, 351.0};

    s_factor = B;
    for (G4int i = 0; i < 4; ++i) {
      s_factor += A[i] / (Square(x_keV - ER[i]) + Square(dE[i]));
    }
  }

  return GetSigmaFromSFactor(energy_cm_MeV, s_factor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BCrossSection::GetSigmaFromSFactor(G4double energy_cm_MeV, G4double s_factor_MeV_b)
{
  if (energy_cm_MeV <= 0. || s_factor_MeV_b <= 0.) return 0.;

  const G4double sigma_b = s_factor_MeV_b / energy_cm_MeV * std::exp(-std::sqrt(GamowEnergyMeV / energy_cm_MeV));
  return sigma_b * barn;
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
