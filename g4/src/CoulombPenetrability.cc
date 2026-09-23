#include "CoulombPenetrability.hh"

#include "CoulombPenetrabilityData.hh"

#include <algorithm>
#include <cmath>

namespace {

// Linear interpolation of log(P_L) against log(E).
// The tables store full Coulomb penetrability P_L(E,a) = rho / (F_L^2 + G_L^2).
double InterpolateLogP(double energy_keV, const double* energy_table, const double* logp_table, int n)
{
  if (energy_keV <= energy_table[0]) return logp_table[0];
  if (energy_keV >= energy_table[n - 1]) return logp_table[n - 1];

  const double* upper = std::upper_bound(energy_table, energy_table + n, energy_keV);
  const int i1 = static_cast<int>(upper - energy_table);
  const int i0 = i1 - 1;

  const double x = std::log(energy_keV);
  const double x0 = std::log(energy_table[i0]);
  const double x1 = std::log(energy_table[i1]);

  const double t = (x - x0) / (x1 - x0);
  return logp_table[i0] + t * (logp_table[i1] - logp_table[i0]);
}

double PenetrabilityRatio(const double energy_keV, const double reference_energy_keV, const double* energy_table, const double* logp_table, const int n)
{
  if (energy_keV <= 0. || reference_energy_keV <= 0.) return 0.;

  const double logp_e = InterpolateLogP(energy_keV, energy_table, logp_table, n);
  const double logp_r = InterpolateLogP(reference_energy_keV, energy_table, logp_table, n);

  return std::exp(logp_e - logp_r);
}

} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace CoulombPenetrability {

double P11BRatio(const double energy_keV, const double reference_energy_keV, const int orbital_l)
{
  if (orbital_l == 0) {
    return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::P11B_energy_keV, CoulombPenetrabilityData::P11B_L0_logP, CoulombPenetrabilityData::P11B_N);
  }

  if (orbital_l == 1) {
    return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::P11B_energy_keV, CoulombPenetrabilityData::P11B_L1_logP, CoulombPenetrabilityData::P11B_N);
  }

  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double AlphaAlphaL2Ratio(const double energy_keV, const double reference_energy_keV)
{
  return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::AlphaAlpha_energy_keV, CoulombPenetrabilityData::AlphaAlpha_L2_logP, CoulombPenetrabilityData::AlphaAlpha_N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Alpha8BeRatio(const double energy_keV, const double reference_energy_keV, const int orbital_l)
{
  if (energy_keV < CoulombPenetrabilityData::Alpha8Be_energy_keV[0]) return 0.;

  if (orbital_l == 1) {
    return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::Alpha8Be_energy_keV, CoulombPenetrabilityData::Alpha8Be_L1_logP, CoulombPenetrabilityData::Alpha8Be_N);
  }

  if (orbital_l == 2) {
    return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::Alpha8Be_energy_keV, CoulombPenetrabilityData::Alpha8Be_L2_logP, CoulombPenetrabilityData::Alpha8Be_N);
  }

  if (orbital_l == 3) {
    return PenetrabilityRatio(energy_keV, reference_energy_keV, CoulombPenetrabilityData::Alpha8Be_energy_keV, CoulombPenetrabilityData::Alpha8Be_L3_logP, CoulombPenetrabilityData::Alpha8Be_N);
  }

  return 0.;
}

} // namespace CoulombPenetrability
