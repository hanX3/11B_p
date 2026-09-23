#ifndef CoulombPenetrability_h
#define CoulombPenetrability_h 1

#include "CoulombPenetrabilityData.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace H11BCoulomb {

// Linear interpolation of log(P_L) against log(E).
// The tables are generated from the full Coulomb penetrability
//     P_L(E,a) = rho / [F_L(eta,rho)^2 + G_L(eta,rho)^2],
// where rho = k a.  Returning ratios avoids changing the partial width at Er.
inline double InterpolateLogP(double energy_keV,
                              const double *energy_table,
                              const double *logp_table,
                              int n)
{
  if(energy_keV <= energy_table[0]) return logp_table[0];
  if(energy_keV >= energy_table[n-1]) return logp_table[n-1];

  const double *upper = std::upper_bound(energy_table, energy_table + n, energy_keV);
  int i1 = static_cast<int>(upper - energy_table);
  int i0 = i1 - 1;

  const double x  = std::log(energy_keV);
  const double x0 = std::log(energy_table[i0]);
  const double x1 = std::log(energy_table[i1]);

  const double t = (x - x0)/(x1 - x0);
  return logp_table[i0] + t*(logp_table[i1] - logp_table[i0]);
}

inline double PenetrabilityRatio(const double energy_keV,
                                 const double reference_energy_keV,
                                 const double *energy_table,
                                 const double *logp_table,
                                 const int n)
{
  if(energy_keV <= 0. || reference_energy_keV <= 0.) return 0.;

  const double logp_e = InterpolateLogP(energy_keV, energy_table, logp_table, n);
  const double logp_r = InterpolateLogP(reference_energy_keV, energy_table, logp_table, n);

  return std::exp(logp_e - logp_r);
}

inline double P11BRatio(const double energy_keV,
                        const double reference_energy_keV,
                        const int orbital_l)
{
  if(orbital_l == 0){
    return PenetrabilityRatio(energy_keV,
                              reference_energy_keV,
                              H11BCoulombData::P11B_energy_keV,
                              H11BCoulombData::P11B_L0_logP,
                              H11BCoulombData::P11B_N);
  }

  if(orbital_l == 1){
    return PenetrabilityRatio(energy_keV,
                              reference_energy_keV,
                              H11BCoulombData::P11B_energy_keV,
                              H11BCoulombData::P11B_L1_logP,
                              H11BCoulombData::P11B_N);
  }

  return 0.;
}

inline double AlphaAlphaL2Ratio(const double energy_keV,
                                const double reference_energy_keV)
{
  return PenetrabilityRatio(energy_keV,
                            reference_energy_keV,
                            H11BCoulombData::AlphaAlpha_energy_keV,
                            H11BCoulombData::AlphaAlpha_L2_logP,
                            H11BCoulombData::AlphaAlpha_N);
}

} // namespace H11BCoulomb

#endif
