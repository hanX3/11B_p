#ifndef CoulombPenetrability_H
#define CoulombPenetrability_H 1

namespace CoulombPenetrability {

double P11BRatio(double energy_keV, double reference_energy_keV, int orbital_l);
double AlphaAlphaL2Ratio(double energy_keV, double reference_energy_keV);
double Alpha8BeRatio(double energy_keV, double reference_energy_keV, int orbital_l);

} // namespace CoulombPenetrability

#endif
