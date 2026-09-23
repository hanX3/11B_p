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
  G4double sigma_165_model = 0.;
  G4double sigma_675_model = 0.;
  G4double sigma_165 = 0.;
  G4double sigma_675 = 0.;
  G4double sigma_background = 0.;
  G4double sigma_165_sampling = 0.;
  G4double sigma_675_sampling = 0.;
  G4double sigma_background_sampling = 0.;
  G4double sigma_3alpha_sampling_total = 0.;
  G4double sigma_eval = 0.;
  G4double sigma_total = 0.;
  G4double sigma_gamma_165_0 = 0.;
  G4double sigma_gamma_165_1 = 0.;
  G4double sigma_gamma_165_total = 0.;
  G4double sigma_gamma_675_total = 0.;
  G4double sigma_gamma_total = 0.;
  G4double sigma_gamma_165_0_physical = 0.;
  G4double sigma_gamma_165_1_physical = 0.;
  G4double sigma_gamma_675_physical = 0.;
  G4double sigma_gamma_total_physical = 0.;
  G4double sigma_gamma_165_0_sampling = 0.;
  G4double sigma_gamma_165_1_sampling = 0.;
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
  G4double background_bias_factor = 1.;
};

enum class H11BEvaluatedCrossSectionMode
{
  Tentori2023,
  Model
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
  static H11BEvaluatedCrossSectionMode GetEvaluatedCrossSectionMode();
  static void SetEvaluatedCrossSectionMode(H11BEvaluatedCrossSectionMode mode);
  static void SetEvaluatedCrossSectionMode(const G4String& mode);
  static const char* GetEvaluatedCrossSectionModeName();
  static G4double GetCrossSectionBiasFactor();
  static void SetCrossSectionBiasFactor(G4double factor);
  static G4double GetBackgroundBiasFactor();
  static void SetBackgroundBiasFactor(G4double factor);

private:
  static G4double GetEcmValue(G4double project_a, G4double target_a, G4double kinetic_energy_lab_keV);
  static G4double GetSigma165(G4double energy_cm_keV);
  static G4double GetSigma675(G4double energy_cm_keV);
  static G4double Get165Gamma0CrossSection(G4double energy_cm_keV);
  static G4double Get165Gamma1CrossSection(G4double energy_cm_keV);
  static G4double Get675GammaTotalCrossSection(G4double sigma_675_model);
  static G4double GetEvaluatedTotalCrossSection(G4double kinetic_energy_lab);
  static G4double GetTentori2023TotalCrossSection(G4double energy_cm_MeV);
  static G4double GetSigmaFromSFactor(G4double energy_cm_MeV, G4double s_factor_MeV_b);
  static G4double GetSigmaBreitWigner(G4double energy_cm_keV, G4double resonance_energy, G4double total_width, G4double entrance_width_at_resonance, G4double exit_width_at_resonance, G4double spin_stat_factor, G4int entrance_orbital_l);
  static G4double GetFullCoulombEntranceWidth(G4double energy_cm_keV, G4double resonance_energy, G4double entrance_width_at_resonance, G4int entrance_orbital_l);
};

#endif
