#ifndef H11BAngularDistribution_h
#define H11BAngularDistribution_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"

namespace H11BAngularDistribution
{
  // Isotropic unit vector in the local frame.
  G4ThreeVector SampleIsotropicDirection();

  // Generic axially-symmetric primary-alpha angular distribution around a
  // chosen symmetry axis.  The weight is
  //   W(theta) = 1 + a2 P2(cos(theta)) + a4 P4(cos(theta)).
  // The returned vector is a unit vector in the same coordinate system as axis.
  G4ThreeVector SampleDirectionFromLegendreA2A4(const G4ThreeVector& axis,
                                                G4double a2,
                                                G4double a4);

  // 165-keV resonance first breakup:
  //   12C*(16.11, 2+) -> alpha + 8Be(g.s./2+).
  // This samples the primary-alpha direction relative to the beam axis in the
  // 12C CM frame.  The alpha0 and alpha1 branches have independent Legendre
  // coefficients in Constants.hh.
  G4ThreeVector Sample165PrimaryAlphaDirection(const G4ThreeVector& beam_axis_cm12c,
                                                G4bool alpha1_branch);

  // 165-keV resonance alpha1 channel:
  //   12C*(16.11, 2+) -> alpha + 8Be(2+), mainly L=2
  //   8Be(2+) -> alpha + alpha, L=2
  // The returned vector is the alpha direction in the 8Be rest frame.
  // It is correlated with the primary alpha direction in the 12C CM frame.
  G4ThreeVector Sample165Alpha1DWaveDirection(const G4ThreeVector& primary_alpha_dir_cm12c);

  // 675-keV resonance alpha1 channel:
  //   12C*(16.57/16.62, 2-) -> alpha + 8Be(2+), coherent L=1/L=3 mixture
  //   8Be(2+) -> alpha + alpha, L=2
  // The returned vector is the alpha direction in the 8Be rest frame.
  // It is correlated with the primary alpha direction in the 12C CM frame.
  G4ThreeVector Sample675Alpha1L13Direction(const G4ThreeVector& primary_alpha_dir_cm12c);

  // Weight functions in x = cos(chi), where chi is the angle between the
  // primary alpha direction and one alpha from the 8Be decay.
  G4double WeightLegendreA2A4(G4double x, G4double a2, G4double a4);
  G4double Weight165Alpha1DWave(G4double x);
  G4double Weight675Alpha1L13(G4double x);
}

#endif
