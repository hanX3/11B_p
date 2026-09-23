#ifndef H11BCrossSection_H
#define H11BCrossSection_H 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct H11BCrossSectionComponents
{
  G4double sigma_162_model = 0.;
  G4double sigma_675_model = 0.;
  G4double sigma_162_total = 0.;
  G4double sigma_675_total = 0.;
  G4double sigma_162 = 0.;
  G4double sigma_675 = 0.;
  G4double sigma_162_directdecay = 0.;
  G4double sigma_675_directdecay = 0.;
  G4double sigma_directdecay = 0.;
  // Backward-compatible alias for old analysis code.
  G4double sigma_background = 0.;
  G4double sigma_162_sampling = 0.;
  G4double sigma_675_sampling = 0.;
  G4double sigma_162_directdecay_sampling = 0.;
  G4double sigma_675_directdecay_sampling = 0.;
  G4double sigma_directdecay_sampling = 0.;
  // Backward-compatible alias for old analysis code.
  G4double sigma_background_sampling = 0.;
  G4double sigma_3alpha_sampling_total = 0.;
  G4double sigma_eval = 0.;
  G4double sigma_total = 0.;
  G4double sigma_gamma_162_0 = 0.;
  G4double sigma_gamma_162_1 = 0.;
  G4double sigma_gamma_162_total = 0.;
  G4double sigma_gamma_675_total = 0.;
  G4double sigma_gamma_total = 0.;
  G4double sigma_gamma_162_0_physical = 0.;
  G4double sigma_gamma_162_1_physical = 0.;
  G4double sigma_gamma_675_physical = 0.;
  G4double sigma_gamma_total_physical = 0.;
  G4double sigma_gamma_162_0_sampling = 0.;
  G4double sigma_gamma_162_1_sampling = 0.;
  G4double sigma_gamma_675_sampling = 0.;
  G4double sigma_gamma_sampling_total = 0.;
  G4double sigma_gamma_675_to_ground_physical = 0.;
  G4double sigma_gamma_675_to_4439_physical = 0.;
  G4double sigma_gamma_675_to_7654_physical = 0.;
  G4double sigma_gamma_675_to_12710_physical = 0.;
  G4double sigma_gamma_675_to_15110_physical = 0.;
  G4double sigma_gamma_675_to_ground_sampling = 0.;
  G4double sigma_gamma_675_to_4439_sampling = 0.;
  G4double sigma_gamma_675_to_7654_sampling = 0.;
  G4double sigma_gamma_675_to_12710_sampling = 0.;
  G4double sigma_gamma_675_to_15110_sampling = 0.;
  G4double sigma_total_physical_all = 0.;
  G4double sigma_total_sampling_all = 0.;
  G4double sigma_total_all = 0.;
  G4double cross_section_bias_factor = 1.;
  G4double direct_decay_fraction = 0.01;
  G4double sequential_decay_fraction_162 = 0.99;
  G4double sequential_decay_fraction_675 = 0.99;
  G4double direct_decay_fraction_162 = 0.01;
  G4double direct_decay_fraction_675 = 0.01;
  G4double scale_factor_162 = 1.;
  G4double scale_factor_675 = 1.;
  G4bool enable_direct_decay = true;
  // Backward-compatible alias for old analysis code.
  G4double background_bias_factor = 0.01;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class H11BCrossSection : public G4VCrossSectionDataSet
{
public:
  H11BCrossSection() = default;
  ~H11BCrossSection() override = default;

  G4bool
  IsIsoApplicable(const G4DynamicParticle* projectile, G4int z, G4int a, const G4Element*, const G4Material*) override;
  G4double GetIsoCrossSection(const G4DynamicParticle* projectile, G4int z, G4int a, const G4Isotope*, const G4Element*, const G4Material*) override;

  static H11BCrossSectionComponents CalculateComponents(G4double kinetic_energy_lab);
  static G4double GetCrossSectionBiasFactor();
  static void SetCrossSectionBiasFactor(G4double factor);
  static G4double Get162SequentialDecayFraction();
  static void Set162SequentialDecayFraction(G4double fraction);
  static G4double Get675SequentialDecayFraction();
  static void Set675SequentialDecayFraction(G4double fraction);
  static G4bool GetDirectDecayEnabled();
  static void SetDirectDecayEnabled(G4bool enabled);
  static G4double Get162BWScaleFactor();
  static void Set162BWScaleFactor(G4double factor);
  static G4double Get675ScaleFactor();
  static void Set675ScaleFactor(G4double factor);

private:
  static G4double GetEcmValue(G4double project_a, G4double target_a, G4double kinetic_energy_lab_keV);
  static G4double GetSigma162(G4double energy_cm_keV);
  static G4double Get162Gamma0CrossSection(G4double energy_cm_keV);
  static G4double Get162Gamma1CrossSection(G4double energy_cm_keV);
  static G4double Get675GammaTotalCrossSection(G4double sigma_675_model);
  static G4double GetFit675CrossSection(G4double kinetic_energy_lab);
  static G4double GetSigmaBreitWigner(G4double energy_cm_keV, G4double resonance_energy, G4double total_width, G4double entrance_width_at_resonance, G4double exit_width_at_resonance, G4double spin_stat_factor, G4int entrance_orbital_l);
  static G4double GetFullCoulombEntranceWidth(G4double energy_cm_keV, G4double resonance_energy, G4double entrance_width_at_resonance, G4int entrance_orbital_l);
};

#endif
