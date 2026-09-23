#ifndef H11BAngularDistribution_H
#define H11BAngularDistribution_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <utility>

namespace H11BAngularDistribution {
// Isotropic unit vector in the local frame.
G4ThreeVector SampleIsotropicDirection();

// Generic axially-symmetric primary-alpha angular distributions around a
// chosen symmetry axis.  The returned vector is a unit vector in the same
// coordinate system as axis.
//   W(theta) = 1 + a1 P1(cos(theta)) + a2 P2(cos(theta))
G4ThreeVector SampleDirectionFromLegendreA1A2(const G4ThreeVector& axis, G4double a1, G4double a2);

// Legacy helper for any future A2/A4 distribution:
//   W(theta) = 1 + a2 P2(cos(theta)) + a4 P4(cos(theta)).
// The returned vector is a unit vector in the same coordinate system as axis.
G4ThreeVector SampleDirectionFromLegendreA2A4(const G4ThreeVector& axis, G4double a2, G4double a4);

// 162-keV resonance first breakup:
//   12C*(16.11, 2+) -> alpha + 8Be(g.s./2+).
// This samples the primary-alpha direction relative to the beam axis in the
// 12C CM frame.  The alpha0 and alpha1 branches have independent Legendre
// coefficients in Constants.hh.
G4ThreeVector
Sample162PrimaryAlphaDirection(const G4ThreeVector& beam_axis_cm12c, G4bool alpha1_branch, G4double ecm_keV);
std::pair<G4double, G4double> Get162PrimaryA1A2(G4double ecm_keV, G4bool alpha1_branch);

// 162-keV resonance alpha1 channel:
//   12C*(16.11, 2+) -> alpha + 8Be(2+), mainly L=2
//   8Be(2+) -> alpha + alpha, L=2
// The returned vector is the alpha direction in the 8Be rest frame.
// It is correlated with the primary alpha direction in the 12C CM frame.
G4ThreeVector Sample162Alpha1DWaveDirection(const G4ThreeVector& primary_alpha_dir_cm12c);

// 675-keV resonance alpha1 channel:
//   12C*(16.57/16.62, 2-) -> alpha + 8Be(2+)
//   8Be(2+) -> alpha + alpha.
// The returned vector is the alpha direction in the 8Be rest frame.
// It is correlated with the 8Be recoil direction in the 12C CM frame.
G4ThreeVector Sample675Alpha1LegendreA2A4Direction(const G4ThreeVector& axis_8be_recoil_cm12c);

// 162-keV ground-state gamma angular distribution in fixed A1/A2 form:
//   W(theta) = 1 + a1*P1(cos(theta)) + a2*P2(cos(theta)),
// where theta is relative to the incident proton beam axis.
G4ThreeVector Sample162Gamma0FixedA1A2Direction(const G4ThreeVector& beam_axis_cm12c);

// Weight functions in x = cos(chi), where chi is the angle between the
// chosen internal axis and one alpha from the 8Be decay.
G4double WeightLegendreA1A2(G4double x, G4double a1, G4double a2);
G4double WeightLegendreA2A4(G4double x, G4double a2, G4double a4);
G4double Weight162Alpha1DWave(G4double x);
G4double Weight675Alpha1LegendreA2A4(G4double x);
G4double Weight162Gamma0FixedA1A2(G4double x);
} // namespace H11BAngularDistribution

#endif
