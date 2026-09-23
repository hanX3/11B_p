#ifndef CoulombPenetrabilityData_H
#define CoulombPenetrabilityData_H 1

// Auto-generated with mpmath Coulomb wave functions.
// The arrays store log(P_L), where P_L = rho/(F_L^2 + G_L^2).
// p + 11B tables: Z1=1, Z2=5, mu=11/12 amu, channel radius a=4.5 fm.
// alpha + alpha table: Z1=2, Z2=2, mu=2 amu, channel radius a=4.5 fm.
// alpha + 8Be table: Z1=2, Z2=4, mu=8/3 amu, channel radius a=4.5 fm.

namespace CoulombPenetrabilityData {

extern const int P11B_N;
extern const int AlphaAlpha_N;
extern const int Alpha8Be_N;

extern const double P11B_energy_keV[];
extern const double P11B_L0_logP[];
extern const double P11B_L1_logP[];

extern const double AlphaAlpha_energy_keV[];
extern const double AlphaAlpha_L2_logP[];

extern const double Alpha8Be_energy_keV[];
extern const double Alpha8Be_L1_logP[];
extern const double Alpha8Be_L2_logP[];
extern const double Alpha8Be_L3_logP[];

} // namespace CoulombPenetrabilityData

#endif
