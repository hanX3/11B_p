#include "H11BAngularDistribution.hh"

#include "Constants.hh"
#include "H11BConfig.hh"

#include "G4Exception.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cmath>
#include <utility>

namespace {
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double ClipNonNegative(G4double value)
{
  static G4bool warned = false;
  if (value < -1e-8 && !warned) {
    G4cerr << "Warning: negative angular weight detected. Check H11B angular-distribution parameters." << G4endl;
    warned = true;
  }

  return value > 0.0 ? value : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector DirectionFromCosChi(const G4ThreeVector& axis, G4double cos_chi)
{
  G4ThreeVector ez = axis;
  if (ez.mag2() <= 0.0) {
    return H11BAngularDistribution::SampleIsotropicDirection();
  }
  ez = ez.unit();

  const G4double sin_chi = std::sqrt(std::max(0.0, 1.0 - cos_chi * cos_chi));
  const G4double phi = twopi * G4UniformRand();

  G4ThreeVector dir_local(sin_chi * std::cos(phi), sin_chi * std::sin(phi), cos_chi);

  // 把局部 z 轴旋转到真正的 axis 方向
  dir_local.rotateUz(ez);

  return dir_local.unit();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double FindWeightMaximum(G4double (*weight_func)(G4double))
{
  // Numeric maximum makes the rejection sampler robust if the user changes
  // the L=1/L=3 fraction or the effective phase in Constants.hh.
  G4double max_weight = 0.0;
  constexpr G4int n_grid = 2000;

  for (G4int i = 0; i <= n_grid; ++i) {
    const G4double x = -1.0 + 2.0 * static_cast<G4double>(i) / n_grid;
    max_weight = std::max(max_weight, ClipNonNegative(weight_func(x)));
  }

  // Avoid zero division and leave a small safety margin.
  return max_weight > 0.0 ? 1.05 * max_weight : 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double SampleCosChi(G4double (*weight_func)(G4double), G4double w_max)
{
  for (;;) {
    const G4double x = 2.0 * G4UniformRand() - 1.0;
    const G4double w = ClipNonNegative(weight_func(x));

    if (G4UniformRand() * w_max <= w) {
      return x;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double FindLegendreWeightMaximum(G4double a2, G4double a4)
{
  G4double max_weight = 0.0;
  constexpr G4int n_grid = 2000;

  for (G4int i = 0; i <= n_grid; ++i) {
    const G4double x = -1.0 + 2.0 * static_cast<G4double>(i) / n_grid;
    max_weight = std::max(max_weight, ClipNonNegative(H11BAngularDistribution::WeightLegendreA2A4(x, a2, a4)));
  }

  return max_weight > 0.0 ? 1.05 * max_weight : 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double SampleCosThetaLegendre(G4double a2, G4double a4, G4double w_max)
{
  for (;;) {
    const G4double x = 2.0 * G4UniformRand() - 1.0;
    const G4double w = ClipNonNegative(H11BAngularDistribution::WeightLegendreA2A4(x, a2, a4));

    if (G4UniformRand() * w_max <= w) {
      return x;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double FindLegendreA1A2Maximum(G4double a1, G4double a2)
{
  G4double max_weight = 0.0;
  constexpr G4int n_grid = 2000;

  for (G4int i = 0; i <= n_grid; ++i) {
    const G4double x = -1.0 + 2.0 * static_cast<G4double>(i) / n_grid;
    const G4double w = H11BAngularDistribution::WeightLegendreA1A2(x, a1, a2);
    if (w < 0.0) {
      G4Exception("H11BAngularDistribution::FindLegendreA1A2Maximum",
                  "H11B162PrimaryAngular001",
                  FatalException,
                  "162-keV primary alpha angular distribution has negative weight.");
    }
    max_weight = std::max(max_weight, w);
  }

  return max_weight > 0.0 ? 1.05 * max_weight : 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double SampleCosThetaLegendreA1A2(G4double a1, G4double a2, G4double w_max)
{
  for (;;) {
    const G4double x = 2.0 * G4UniformRand() - 1.0;
    const G4double w = H11BAngularDistribution::WeightLegendreA1A2(x, a1, a2);
    if (w < 0.0) {
      G4Exception("H11BAngularDistribution::SampleCosThetaLegendreA1A2",
                  "H11B162PrimaryAngular001",
                  FatalException,
                  "162-keV primary alpha angular distribution has negative weight.");
    }

    if (G4UniformRand() * w_max <= w) {
      return x;
    }
  }
}

} // namespace

namespace H11BAngularDistribution {
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector SampleIsotropicDirection()
{
  const G4double cos_theta = 2.0 * G4UniformRand() - 1.0;
  const G4double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
  const G4double phi = twopi * G4UniformRand();

  return G4ThreeVector(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double WeightLegendreA1A2(G4double x, G4double a1, G4double a2)
{
  const G4double p1 = x;
  const G4double p2 = 0.5 * (3.0 * x * x - 1.0);

  return 1.0 + a1 * p1 + a2 * p2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double WeightLegendreA2A4(G4double x, G4double a2, G4double a4)
{
  const G4double x2 = x * x;
  const G4double x4 = x2 * x2;

  const G4double p2 = 0.5 * (3.0 * x2 - 1.0);
  const G4double p4 = (1.0 / 8.0) * (35.0 * x4 - 30.0 * x2 + 3.0);

  return 1.0 + a2 * p2 + a4 * p4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector SampleDirectionFromLegendreA1A2(const G4ThreeVector& axis, G4double a1, G4double a2)
{
  if (axis.mag2() <= 0.0) {
    return SampleIsotropicDirection();
  }

  const G4double w_max = FindLegendreA1A2Maximum(a1, a2);
  const G4double cos_theta = SampleCosThetaLegendreA1A2(a1, a2, w_max);

  return DirectionFromCosChi(axis, cos_theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector SampleDirectionFromLegendreA2A4(const G4ThreeVector& axis, G4double a2, G4double a4)
{
  if (axis.mag2() <= 0.0) {
    return SampleIsotropicDirection();
  }

  const G4double w_max = FindLegendreWeightMaximum(a2, a4);
  const G4double cos_theta = SampleCosThetaLegendre(a2, a4, w_max);

  return DirectionFromCosChi(axis, cos_theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::pair<G4double, G4double> Get162PrimaryA1A2(G4double ecm_keV, G4bool alpha1_branch)
{
  (void)ecm_keV;

  if (!H11BConfig::GetEnable162PrimaryAngularDistribution()) {
    return {0.0, 0.0};
  }

  return H11BConfig::Get162PrimaryA1A2(alpha1_branch);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector
Sample162PrimaryAlphaDirection(const G4ThreeVector& beam_axis_cm12c, G4bool alpha1_branch, G4double ecm_keV)
{
  if (beam_axis_cm12c.mag2() <= 0.0 ||
      !H11BConfig::GetEnable162PrimaryAngularDistribution()) {
    return SampleIsotropicDirection();
  }

  const auto [a1, a2] = Get162PrimaryA1A2(ecm_keV, alpha1_branch);

  return SampleDirectionFromLegendreA1A2(beam_axis_cm12c, a1, a2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double Weight162Alpha1DWave(G4double x)
{
  const G4double x2 = x * x;
  const G4double x4 = x2 * x2;

  // Sequential angular correlation for an un-oriented 2+ parent:
  //   2+ -> 2+ + alpha with L=2,
  //   2+ -> 0+ + 0+ with L=2.
  // Normalized over x in [-1,1].
  return (5.0 / 28.0) * (9.0 * x4 - 9.0 * x2 + 4.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double Weight675Alpha1L13(G4double x)
{
  const G4double x2 = x * x;
  const G4double x4 = x2 * x2;

  // Pure L=1 and L=3 components for the 2- -> 2+ -> 0+ cascade.
  const G4double P1 = 0.75 * (1.0 - x2);
  const G4double P3 = (3.0 / 16.0) * (1.0 + 14.0 * x2 - 15.0 * x4);

  // L=1/L=3 angular interference basis.  Its integral over x is zero.
  const G4double C13 = (3.0 / 8.0) * (1.0 - 6.0 * x2 + 5.0 * x4);

  const G4double k = std::clamp(H11B675ExitL1Fraction, 0.0, 1.0);

  G4double weight = k * P1 + (1.0 - k) * P3;

  if (H11B675UseCoherentL13Interference) {
    weight += 2.0 * std::sqrt(std::max(0.0, k * (1.0 - k))) * std::cos(H11B675ExitL13Phase) * C13;
  }

  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double Weight675Alpha1LegendreA2A4(G4double x)
{
  return WeightLegendreA2A4(x,
                            H11BConfig::Get675Alpha1SecondaryA2(),
                            H11BConfig::Get675Alpha1SecondaryA4());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double Weight162Gamma0FixedA1A2(G4double x)
{
  return WeightLegendreA1A2(x, H11BConfig::Get162Gamma0AngularA1(), H11BConfig::Get162Gamma0AngularA2());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector Sample162Alpha1DWaveDirection(const G4ThreeVector& primary_alpha_dir_cm12c)
{
  static const G4double w_max = FindWeightMaximum(Weight162Alpha1DWave);
  const G4double cos_chi = SampleCosChi(Weight162Alpha1DWave, w_max);
  return DirectionFromCosChi(primary_alpha_dir_cm12c, cos_chi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector Sample675Alpha1LegendreA2A4Direction(const G4ThreeVector& axis_8be_recoil_cm12c)
{
  const G4double w_max = FindWeightMaximum(Weight675Alpha1LegendreA2A4);
  const G4double cos_chi = SampleCosChi(Weight675Alpha1LegendreA2A4, w_max);
  return DirectionFromCosChi(axis_8be_recoil_cm12c, cos_chi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector Sample162Gamma0FixedA1A2Direction(const G4ThreeVector& beam_axis_cm12c)
{
  const G4double w_max = FindWeightMaximum(Weight162Gamma0FixedA1A2);
  const G4double cos_theta = SampleCosChi(Weight162Gamma0FixedA1A2, w_max);
  return DirectionFromCosChi(beam_axis_cm12c, cos_theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreeVector Sample675Alpha1L13Direction(const G4ThreeVector& primary_alpha_dir_cm12c)
{
  static const G4double w_max = FindWeightMaximum(Weight675Alpha1L13);
  const G4double cos_chi = SampleCosChi(Weight675Alpha1L13, w_max);
  return DirectionFromCosChi(primary_alpha_dir_cm12c, cos_chi);
}
} // namespace H11BAngularDistribution
