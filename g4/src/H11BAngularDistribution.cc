#include "H11BAngularDistribution.hh"

#include "Constants.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cmath>

namespace
{
  G4double ClipNonNegative(G4double value)
  {
    return value > 0.0 ? value : 0.0;
  }

  G4ThreeVector DirectionFromCosChi(const G4ThreeVector& axis, G4double cos_chi)
  {
    G4ThreeVector ez = axis;
    if(ez.mag2() <= 0.0){
      return H11BAngularDistribution::SampleIsotropicDirection();
    }
    ez = ez.unit();

    const G4double sin_chi = std::sqrt(std::max(0.0, 1.0 - cos_chi*cos_chi));
    const G4double phi = twopi * G4UniformRand();

    G4ThreeVector dir_local(
      sin_chi * std::cos(phi),
      sin_chi * std::sin(phi),
      cos_chi
  );

  // 把局部 z 轴旋转到真正的 axis 方向
  dir_local.rotateUz(ez);

    return dir_local.unit();
  }

  G4double FindWeightMaximum(G4double (*weight_func)(G4double))
  {
    // Numeric maximum makes the rejection sampler robust if the user changes
    // the L=1/L=3 fraction or the effective phase in Constants.hh.
    G4double max_weight = 0.0;
    constexpr G4int n_grid = 2000;

    for(G4int i=0; i<=n_grid; ++i){
      const G4double x = -1.0 + 2.0 * static_cast<G4double>(i) / n_grid;
      max_weight = std::max(max_weight, ClipNonNegative(weight_func(x)));
    }

    // Avoid zero division and leave a small safety margin.
    return max_weight > 0.0 ? 1.05 * max_weight : 1.0;
  }

  G4double SampleCosChi(G4double (*weight_func)(G4double), G4double w_max)
  {
    for(;;){
      const G4double x = 2.0 * G4UniformRand() - 1.0;
      const G4double w = ClipNonNegative(weight_func(x));

      if(G4UniformRand() * w_max <= w){
        return x;
      }
    }
  }

  G4double FindLegendreWeightMaximum(G4double a2, G4double a4)
  {
    G4double max_weight = 0.0;
    constexpr G4int n_grid = 2000;

    for(G4int i=0; i<=n_grid; ++i){
      const G4double x = -1.0 + 2.0 * static_cast<G4double>(i) / n_grid;
      max_weight = std::max(max_weight,
                            ClipNonNegative(H11BAngularDistribution::WeightLegendreA2A4(x, a2, a4)));
    }

    return max_weight > 0.0 ? 1.05 * max_weight : 1.0;
  }

  G4double SampleCosThetaLegendre(G4double a2, G4double a4, G4double w_max)
  {
    for(;;){
      const G4double x = 2.0 * G4UniformRand() - 1.0;
      const G4double w = ClipNonNegative(H11BAngularDistribution::WeightLegendreA2A4(x, a2, a4));

      if(G4UniformRand() * w_max <= w){
        return x;
      }
    }
  }
}

namespace H11BAngularDistribution
{
  G4ThreeVector SampleIsotropicDirection()
  {
    const G4double cos_theta = 2.0 * G4UniformRand() - 1.0;
    const G4double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta*cos_theta));
    const G4double phi = twopi * G4UniformRand();

    return G4ThreeVector(sin_theta * std::cos(phi),
                         sin_theta * std::sin(phi),
                         cos_theta);
  }

  G4double WeightLegendreA2A4(G4double x, G4double a2, G4double a4)
  {
    const G4double x2 = x*x;
    const G4double x4 = x2*x2;

    const G4double p2 = 0.5 * (3.0*x2 - 1.0);
    const G4double p4 = (1.0/8.0) * (35.0*x4 - 30.0*x2 + 3.0);

    return 1.0 + a2*p2 + a4*p4;
  }

  G4ThreeVector SampleDirectionFromLegendreA2A4(const G4ThreeVector& axis,
                                                G4double a2,
                                                G4double a4)
  {
    if(axis.mag2() <= 0.0){
      return SampleIsotropicDirection();
    }

    const G4double w_max = FindLegendreWeightMaximum(a2, a4);
    const G4double cos_theta = SampleCosThetaLegendre(a2, a4, w_max);

    return DirectionFromCosChi(axis, cos_theta);
  }

  G4ThreeVector Sample165PrimaryAlphaDirection(const G4ThreeVector& beam_axis_cm12c,
                                                G4bool alpha1_branch)
  {
    if(!H11B165UsePrimaryAngularDistribution || beam_axis_cm12c.mag2() <= 0.0){
      return SampleIsotropicDirection();
    }

    const G4double a2 = alpha1_branch ? H11B165Alpha1PrimaryA2 : H11B165Alpha0PrimaryA2;
    const G4double a4 = alpha1_branch ? H11B165Alpha1PrimaryA4 : H11B165Alpha0PrimaryA4;

    // Reuse the generic axially-symmetric Legendre sampler here, so the 165-keV
    // primary-alpha distribution is handled by the same tested path as any
    // future W(theta)=1+a2*P2+a4*P4 distribution.
    return SampleDirectionFromLegendreA2A4(beam_axis_cm12c, a2, a4);
  }

  G4double Weight165Alpha1DWave(G4double x)
  {
    const G4double x2 = x*x;
    const G4double x4 = x2*x2;

    // Sequential angular correlation for an un-oriented 2+ parent:
    //   2+ -> 2+ + alpha with L=2,
    //   2+ -> 0+ + 0+ with L=2.
    // Normalized over x in [-1,1].
    return (5.0/28.0) * (9.0*x4 - 9.0*x2 + 4.0);
  }

  G4double Weight675Alpha1L13(G4double x)
  {
    const G4double x2 = x*x;
    const G4double x4 = x2*x2;

    // Pure L=1 and L=3 components for the 2- -> 2+ -> 0+ cascade.
    const G4double P1 = 0.75 * (1.0 - x2);
    const G4double P3 = (3.0/16.0) * (1.0 + 14.0*x2 - 15.0*x4);

    // L=1/L=3 angular interference basis.  Its integral over x is zero.
    const G4double C13 = (3.0/8.0) * (1.0 - 6.0*x2 + 5.0*x4);

    const G4double k = std::clamp(H11B675ExitL1Fraction, 0.0, 1.0);

    G4double weight = k*P1 + (1.0-k)*P3;

    if(H11B675UseCoherentL13Interference){
      weight += 2.0 * std::sqrt(std::max(0.0, k*(1.0-k)))
              * std::cos(H11B675ExitL13Phase)
              * C13;
    }

    return weight;
  }

  G4ThreeVector Sample165Alpha1DWaveDirection(const G4ThreeVector& primary_alpha_dir_cm12c)
  {
    static const G4double w_max = FindWeightMaximum(Weight165Alpha1DWave);
    const G4double cos_chi = SampleCosChi(Weight165Alpha1DWave, w_max);
    return DirectionFromCosChi(primary_alpha_dir_cm12c, cos_chi);
  }

  G4ThreeVector Sample675Alpha1L13Direction(const G4ThreeVector& primary_alpha_dir_cm12c)
  {
    static const G4double w_max = FindWeightMaximum(Weight675Alpha1L13);
    const G4double cos_chi = SampleCosChi(Weight675Alpha1L13, w_max);
    return DirectionFromCosChi(primary_alpha_dir_cm12c, cos_chi);
  }
}
