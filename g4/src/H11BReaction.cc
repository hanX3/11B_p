#include "H11BReaction.hh"
#include "H11BAngularDistribution.hh"
#include "H11BConfig.hh"
#include "H11BCrossSection.hh"

#include "globals.hh"
#include "G4RunManager.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4GenericMessenger.hh"
#include "G4Gamma.hh"
#include "G4Threading.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <map>
#include <limits>
#include <utility>
#include <vector>

#include "Constants.hh"
#include "CoulombPenetrability.hh"

#include "Randomize.hh"

namespace {
struct ExcitationCdfTable
{
  G4double ex_max = 0.;
  std::vector<G4double> x;
  std::vector<G4double> cdf;
};

struct H11BGammaLine
{
  H11BGammaBranch branch;
  G4double final_state_energy;
  G4double relative_intensity;
  G4double relative_intensity_uncertainty;
  G4bool is_upper_limit;
  const char* label;
};

constexpr H11BGammaLine kH11B675GammaRelativeLines[] = {{H11BGammaBranch::Gamma675ToGroundState, C12LevelGround, 15.7, 1.6, false, "toGround"}, {H11BGammaBranch::Gamma675To4439State, C12Level4439, 100.0, 0.0, false, "to4439"}, {H11BGammaBranch::Gamma675To7654State, C12Level7654, 0.07, 0.0, true, "to7654UpperLimit"}, {H11BGammaBranch::Gamma675To12710State, C12Level12710, 6.8, 0.4, false, "to12710"}, {H11BGammaBranch::Gamma675To15110State, C12Level15110, 0.16, 0.03, false, "to15110"}};

G4double SafeSqrt(G4double x)
{
  return std::sqrt(std::max(0.0, x));
}

G4double SafeCosBetween(const G4ThreeVector& a, const G4ThreeVector& b)
{
  return a.mag2() > 0.0 && b.mag2() > 0.0 ? a.unit().dot(b.unit()) : 0.0;
}

G4double SafeOpeningAngle(const G4ThreeVector& a, const G4ThreeVector& b)
{
  return a.mag2() > 0.0 && b.mag2() > 0.0 ? a.angle(b) : 0.0;
}

G4double AzimuthAroundAxis(const G4ThreeVector& direction, const G4ThreeVector& axis)
{
  if (direction.mag2() <= 0.0 || axis.mag2() <= 0.0) {
    return 0.0;
  }

  const G4ThreeVector ez = axis.unit();
  const G4ThreeVector reference = std::abs(ez.z()) < 0.95 ? G4ThreeVector(0.0, 0.0, 1.0) : G4ThreeVector(0.0, 1.0, 0.0);
  const G4ThreeVector ex = reference.cross(ez).unit();
  const G4ThreeVector ey = ez.cross(ex).unit();
  const G4ThreeVector transverse = direction.unit() - direction.unit().dot(ez) * ez;

  if (transverse.mag2() <= 0.0) {
    return 0.0;
  }

  G4double phi = std::atan2(transverse.dot(ey), transverse.dot(ex));
  if (phi < 0.0) {
    phi += twopi;
  }

  return phi;
}

const H11BGammaLine* Find675GammaLine(H11BGammaBranch branch)
{
  for (const auto& line : kH11B675GammaRelativeLines) {
    if (line.branch == branch) {
      return &line;
    }
  }

  return nullptr;
}

G4double Sum675GammaRelativeIntensity()
{
  G4double sum = 0.0;
  for (const auto& line : kH11B675GammaRelativeLines) {
    if (line.is_upper_limit && !H11B675IncludeUpperLimitGammaLines) {
      continue;
    }
    sum += line.relative_intensity;
  }

  return sum;
}

H11BGammaBranch Sample675GammaRelativeLine()
{
  const G4double sum = Sum675GammaRelativeIntensity();
  if (sum <= 0.0) {
    return H11BGammaBranch::Gamma675To4439State;
  }

  G4double r = G4UniformRand() * sum;
  for (const auto& line : kH11B675GammaRelativeLines) {
    if (line.is_upper_limit && !H11B675IncludeUpperLimitGammaLines) {
      continue;
    }

    r -= line.relative_intensity;
    if (r <= 0.0) {
      return line.branch;
    }
  }

  return H11BGammaBranch::Gamma675To15110State;
}

G4double GammaBranchPhysicalWeight(H11BGammaBranch branch, const H11BCrossSectionComponents& components)
{
  switch (branch) {
  case H11BGammaBranch::Gamma165ToGroundState:
    return components.sigma_gamma_165_0_physical;
  case H11BGammaBranch::Gamma165ToFirstExcitedState:
    return components.sigma_gamma_165_1_physical;
  case H11BGammaBranch::Gamma675ToGroundState:
    return components.sigma_gamma_675_to_ground_physical;
  case H11BGammaBranch::Gamma675To4439State:
    return components.sigma_gamma_675_to_4439_physical;
  case H11BGammaBranch::Gamma675To7654State:
    return components.sigma_gamma_675_to_7654_physical;
  case H11BGammaBranch::Gamma675To12710State:
    return components.sigma_gamma_675_to_12710_physical;
  case H11BGammaBranch::Gamma675To15110State:
    return components.sigma_gamma_675_to_15110_physical;
  case H11BGammaBranch::None:
    return 0.0;
  }

  return 0.0;
}

void FillCrossSectionDiagnostics(H11BReactionData& data, const H11BCrossSectionComponents& components)
{
  data.sigma_eval_b = components.sigma_eval / barn;
  data.sigma_165_model_b = components.sigma_165_model / barn;
  data.sigma_675_model_b = components.sigma_675_model / barn;
  data.sigma_165_used_b = components.sigma_165 / barn;
  data.sigma_675_used_b = components.sigma_675 / barn;
  data.sigma_total_used_b = components.sigma_total / barn;
  data.sigma_165_sampling_b = components.sigma_165_sampling / barn;
  data.sigma_675_sampling_b = components.sigma_675_sampling / barn;
  data.sigma_background_sampling_b = components.sigma_background_sampling / barn;
  data.sigma_3alpha_sampling_total_b = components.sigma_3alpha_sampling_total / barn;
  data.sigma_model_sum_b = (components.sigma_165_model + components.sigma_675_model) / barn;
  data.model_scale_factor = data.sigma_model_sum_b > 0.0 ? (components.sigma_165 + components.sigma_675) / (components.sigma_165_model + components.sigma_675_model) : 1.0;
  data.cross_section_bias_factor = components.cross_section_bias_factor;
  data.background_bias_factor = components.background_bias_factor;
  data.sigma_background_b = components.sigma_background / barn;
  data.sigma_3alpha_eval_b = components.sigma_eval / barn;
  data.sigma_gamma_165_0_b = components.sigma_gamma_165_0 / barn;
  data.sigma_gamma_165_1_b = components.sigma_gamma_165_1 / barn;
  data.sigma_gamma_165_total_b = components.sigma_gamma_165_total / barn;
  data.sigma_gamma_675_total_b = components.sigma_gamma_675_total / barn;
  data.sigma_gamma_total_b = components.sigma_gamma_total / barn;
  data.sigma_total_physical_all_b = components.sigma_total_physical_all / barn;
  data.sigma_total_sampling_all_b = components.sigma_total_sampling_all / barn;
  data.sigma_gamma_165_0_physical_b = components.sigma_gamma_165_0_physical / barn;
  data.sigma_gamma_165_1_physical_b = components.sigma_gamma_165_1_physical / barn;
  data.sigma_gamma_675_physical_b = components.sigma_gamma_675_physical / barn;
  data.sigma_gamma_165_0_sampling_b = components.sigma_gamma_165_0_sampling / barn;
  data.sigma_gamma_165_1_sampling_b = components.sigma_gamma_165_1_sampling / barn;
  data.sigma_gamma_675_sampling_b = components.sigma_gamma_675_sampling / barn;
  data.sigma_gamma_675_to_ground_physical_b = components.sigma_gamma_675_to_ground_physical / barn;
  data.sigma_gamma_675_to_4439_physical_b = components.sigma_gamma_675_to_4439_physical / barn;
  data.sigma_gamma_675_to_7654_physical_b = components.sigma_gamma_675_to_7654_physical / barn;
  data.sigma_gamma_675_to_12710_physical_b = components.sigma_gamma_675_to_12710_physical / barn;
  data.sigma_gamma_675_to_15110_physical_b = components.sigma_gamma_675_to_15110_physical / barn;
  data.sigma_gamma_675_to_ground_sampling_b = components.sigma_gamma_675_to_ground_sampling / barn;
  data.sigma_gamma_675_to_4439_sampling_b = components.sigma_gamma_675_to_4439_sampling / barn;
  data.sigma_gamma_675_to_7654_sampling_b = components.sigma_gamma_675_to_7654_sampling / barn;
  data.sigma_gamma_675_to_12710_sampling_b = components.sigma_gamma_675_to_12710_sampling / barn;
  data.sigma_gamma_675_to_15110_sampling_b = components.sigma_gamma_675_to_15110_sampling / barn;
  data.sigma_gamma_sampling_total_b = components.sigma_gamma_sampling_total / barn;
  data.sigma_total_all_b = components.sigma_total_all / barn;

  if (components.sigma_total > 0.) {
    data.channel_probability_165 = components.sigma_165 / components.sigma_total;
    data.channel_probability_675 = components.sigma_675 / components.sigma_total;
    data.channel_probability_background = components.sigma_background / components.sigma_total;
  }

  if (components.sigma_total_sampling_all > 0.) {
    data.channel_probability_gamma = components.sigma_gamma_sampling_total / components.sigma_total_sampling_all;
    data.probability_gamma_165_0 = components.sigma_gamma_165_0_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_165_1 = components.sigma_gamma_165_1_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_total = components.sigma_gamma_675_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_to_ground = components.sigma_gamma_675_to_ground_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_to_4439 = components.sigma_gamma_675_to_4439_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_to_7654 = components.sigma_gamma_675_to_7654_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_to_12710 = components.sigma_gamma_675_to_12710_sampling / components.sigma_total_sampling_all;
    data.probability_gamma_675_to_15110 = components.sigma_gamma_675_to_15110_sampling / components.sigma_total_sampling_all;
  }
}

void FillRuntimeConfigDiagnostics(H11BReactionData& data)
{
  data.enable_165_primary_angular_distribution = H11BConfig::GetEnable165PrimaryAngularDistribution() ? 1 : 0;
  data.a1_165_primary = H11BConfig::Get165PrimaryAngularA1();
  data.a2_165_primary = H11BConfig::Get165PrimaryAngularA2();
  data.enable_675_primary_angular_distribution = H11BConfig::GetEnable675PrimaryAngularDistribution() ? 1 : 0;
  data.a1_675_primary = H11BConfig::Get675PrimaryAngularA1();
  data.a2_675_primary = H11BConfig::Get675PrimaryAngularA2();
  data.enable_675_alpha1_secondary_angular_correlation =
    H11BConfig::GetEnable675Alpha1SecondaryAngularCorrelation() ? 1 : 0;
  data.enable_675_gamma_angular_distribution = H11BConfig::GetEnable675GammaAngularDistribution() ? 1 : 0;
  data.a1_675_gamma = H11BConfig::Get675GammaAngularA1();
  data.a2_675_gamma = H11BConfig::Get675GammaAngularA2();
}

G4double ReactionChannelPhysicalWeight(H11BReaction::ReactionChannel channel, const H11BCrossSectionComponents& components)
{
  switch (channel) {
  case H11BReaction::ReactionChannel::Resonance165:
    return components.sigma_165;
  case H11BReaction::ReactionChannel::Resonance675:
    return components.sigma_675;
  case H11BReaction::ReactionChannel::Background3Alpha:
    return components.sigma_background;
  case H11BReaction::ReactionChannel::GammaCapture12C:
    return components.sigma_gamma_sampling_total > 0.0 ? components.sigma_gamma_total_physical : 0.0;
  }

  return 0.0;
}

G4double ReactionChannelSamplingWeight(H11BReaction::ReactionChannel channel, const H11BCrossSectionComponents& components)
{
  switch (channel) {
  case H11BReaction::ReactionChannel::Resonance165:
    return components.sigma_165_sampling;
  case H11BReaction::ReactionChannel::Resonance675:
    return components.sigma_675_sampling;
  case H11BReaction::ReactionChannel::Background3Alpha:
    return components.sigma_background_sampling;
  case H11BReaction::ReactionChannel::GammaCapture12C:
    return components.sigma_gamma_sampling_total;
  }

  return 0.0;
}

G4double ReactionChannelEventWeight(H11BReaction::ReactionChannel channel, const H11BCrossSectionComponents& components)
{
  const G4double physical = ReactionChannelPhysicalWeight(channel, components);
  const G4double sampling = ReactionChannelSamplingWeight(channel, components);
  return sampling > 0.0 ? physical / sampling : 1.0;
}

G4double GammaBranchWeight(H11BGammaBranch branch, const H11BCrossSectionComponents& components)
{
  switch (branch) {
  case H11BGammaBranch::Gamma165ToGroundState:
    return components.sigma_gamma_165_0_sampling;
  case H11BGammaBranch::Gamma165ToFirstExcitedState:
    return components.sigma_gamma_165_1_sampling;
  case H11BGammaBranch::Gamma675ToGroundState:
    return components.sigma_gamma_675_to_ground_sampling;
  case H11BGammaBranch::Gamma675To4439State:
    return components.sigma_gamma_675_to_4439_sampling;
  case H11BGammaBranch::Gamma675To7654State:
    return components.sigma_gamma_675_to_7654_sampling;
  case H11BGammaBranch::Gamma675To12710State:
    return components.sigma_gamma_675_to_12710_sampling;
  case H11BGammaBranch::Gamma675To15110State:
    return components.sigma_gamma_675_to_15110_sampling;
  case H11BGammaBranch::None:
    return 0.0;
  }

  return 0.0;
}

H11BReaction::BackgroundMode g_background_mode = H11BReaction::BackgroundMode::PhaseSpace;

H11BGammaBranch Select675GammaBranch()
{
  return Sample675GammaRelativeLine();
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H11BReaction::H11BReaction()
    : G4HadronicInteraction()
{
  SetMinEnergy(0. * CLHEP::keV);
  SetMaxEnergy(100. * CLHEP::MeV);
  isBlocked = false;
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H11BReaction::~H11BReaction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
H11BReaction::BackgroundMode H11BReaction::GetBackgroundMode()
{
  return g_background_mode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetBackgroundMode(BackgroundMode mode)
{
  g_background_mode = mode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetBackgroundMode(const G4String& mode)
{
  if (mode == "phaseSpace") {
    SetBackgroundMode(BackgroundMode::PhaseSpace);
  } else {
    G4cerr << "Unknown /h11b/backgroundMode '" << mode << "'. Current version only supports phaseSpace." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const char* H11BReaction::GetBackgroundModeName()
{
  switch (g_background_mode) {
  case BackgroundMode::PhaseSpace:
    return "phaseSpace";
  }

  return "unknown";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEvaluatedCrossSectionModeCommand(const G4String& mode)
{
  H11BCrossSection::SetEvaluatedCrossSectionMode(mode);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetCrossSectionBiasFactorCommand(G4double factor)
{
  H11BCrossSection::SetCrossSectionBiasFactor(factor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetBackgroundBiasFactorCommand(G4double factor)
{
  H11BCrossSection::SetBackgroundBiasFactor(factor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetBackgroundModeCommand(const G4String& mode)
{
  SetBackgroundMode(mode);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable675PrimaryAngularDistributionCommand(G4bool enabled)
{
  H11BConfig::SetEnable675PrimaryAngularDistribution(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675PrimaryAngularA1Command(G4double value)
{
  H11BConfig::Set675PrimaryAngularA1(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675PrimaryAngularA2Command(G4double value)
{
  H11BConfig::Set675PrimaryAngularA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable675Alpha1SecondaryAngularCorrelationCommand(G4bool enabled)
{
  H11BConfig::SetEnable675Alpha1SecondaryAngularCorrelation(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675Alpha1SecondaryA2Command(G4double value)
{
  H11BConfig::Set675Alpha1SecondaryA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675Alpha1SecondaryA4Command(G4double value)
{
  H11BConfig::Set675Alpha1SecondaryA4(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable165PrimaryAngularDistributionCommand(G4bool enabled)
{
  H11BConfig::SetEnable165PrimaryAngularDistribution(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165PrimaryAngularA1Command(G4double value)
{
  H11BConfig::Set165PrimaryAngularA1(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165PrimaryAngularA2Command(G4double value)
{
  H11BConfig::Set165PrimaryAngularA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable165Alpha1SecondaryAngularCorrelationCommand(G4bool enabled)
{
  H11BConfig::SetEnable165Alpha1SecondaryAngularCorrelation(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165Alpha1SecondaryA2Command(G4double value)
{
  H11BConfig::Set165Alpha1SecondaryA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165Alpha1SecondaryA4Command(G4double value)
{
  H11BConfig::Set165Alpha1SecondaryA4(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable165Gamma0AngularDistributionCommand(G4bool enabled)
{
  H11BConfig::SetEnable165Gamma0AngularDistribution(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165Gamma0AngularA1Command(G4double value)
{
  H11BConfig::Set165Gamma0AngularA1(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165Gamma0AngularA2Command(G4double value)
{
  H11BConfig::Set165Gamma0AngularA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetEnable675GammaAngularDistributionCommand(G4bool enabled)
{
  H11BConfig::SetEnable675GammaAngularDistribution(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675GammaAngularA1Command(G4double value)
{
  H11BConfig::Set675GammaAngularA1(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675GammaAngularA2Command(G4double value)
{
  H11BConfig::Set675GammaAngularA2(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set165GammaCaptureEnabledCommand(G4bool enabled)
{
  H11BConfig::Set165GammaCaptureEnabled(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Set675GammaCaptureEnabledCommand(G4bool enabled)
{
  H11BConfig::Set675GammaCaptureEnabled(enabled);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::SetGammaBiasFactorCommand(G4double factor)
{
  H11BConfig::SetGammaBiasFactor(factor);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::PrintConfigCommand()
{
  if (!G4Threading::IsMasterThread()) {
    return;
  }

  G4cout << "H11B reaction configuration" << G4endl
         << "  evaluatedCrossSection        = " << H11BCrossSection::GetEvaluatedCrossSectionModeName() << G4endl
         << "  crossSectionBiasFactor       = " << H11BCrossSection::GetCrossSectionBiasFactor() << G4endl
         << "  backgroundBiasFactor         = " << H11BCrossSection::GetBackgroundBiasFactor() << G4endl
         << "  backgroundMode               = " << GetBackgroundModeName() << G4endl
         << G4endl
         << "  165 primary alpha angular distribution" << G4endl
         << "    enable = " << (H11BConfig::GetEnable165PrimaryAngularDistribution() ? "true" : "false") << G4endl
         << "    a1     = " << H11BConfig::Get165PrimaryAngularA1() << G4endl
         << "    a2     = " << H11BConfig::Get165PrimaryAngularA2() << G4endl
         << "  165 alpha1 secondary angular correlation" << G4endl
         << "    enable = " << (H11BConfig::GetEnable165Alpha1SecondaryAngularCorrelation() ? "true" : "false") << G4endl
         << "    a2     = " << H11BConfig::Get165Alpha1SecondaryA2() << G4endl
         << "    a4     = " << H11BConfig::Get165Alpha1SecondaryA4() << G4endl
         << "  165 gamma0 angular distribution" << G4endl
         << "    enable = " << (H11BConfig::GetEnable165Gamma0AngularDistribution() ? "true" : "false") << G4endl
         << "    a1     = " << H11BConfig::Get165Gamma0AngularA1() << G4endl
         << "    a2     = " << H11BConfig::Get165Gamma0AngularA2() << G4endl
         << G4endl
         << "  675 primary alpha angular distribution" << G4endl
         << "    enable = " << (H11BConfig::GetEnable675PrimaryAngularDistribution() ? "true" : "false") << G4endl
         << "    a1     = " << H11BConfig::Get675PrimaryAngularA1() << G4endl
         << "    a2     = " << H11BConfig::Get675PrimaryAngularA2() << G4endl
         << "  675 alpha1 secondary angular correlation" << G4endl
         << "    enable = " << (H11BConfig::GetEnable675Alpha1SecondaryAngularCorrelation() ? "true" : "false") << G4endl
         << "    a2     = " << H11BConfig::Get675Alpha1SecondaryA2() << G4endl
         << "    a4     = " << H11BConfig::Get675Alpha1SecondaryA4() << G4endl
         << "  675Alpha0Width               = 0 keV (parity forbidden)" << G4endl
         << G4endl
         << "  enable165GammaCapture        = " << (H11BConfig::Get165GammaCaptureEnabled() ? "true" : "false") << G4endl
         << "  enable675GammaCapture        = " << (H11BConfig::Get675GammaCaptureEnabled() ? "true" : "false") << G4endl
         << "  675 gamma angular distribution" << G4endl
         << "    enable = " << (H11BConfig::GetEnable675GammaAngularDistribution() ? "true" : "false") << G4endl
         << "    a1     = " << H11BConfig::Get675GammaAngularA1() << G4endl
         << "    a2     = " << H11BConfig::Get675GammaAngularA2() << G4endl
         << "  675GammaBranchSelection      = relativeLineTable (fixed)" << G4endl
         << "  gammaBiasFactor              = " << H11BConfig::GetGammaBiasFactor() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/h11b/", "p + 11B reaction controls");

  auto& background_mode_cmd = messenger->DeclareMethod("backgroundMode", &H11BReaction::SetBackgroundModeCommand, "Set background 3-alpha mode");
  background_mode_cmd.SetCandidates("phaseSpace");

  messenger->DeclareMethod("backgroundBiasFactor", &H11BReaction::SetBackgroundBiasFactorCommand,
                           "Set background phase-space sampling bias factor (>= 0; 1 disables biasing)");

  auto& xs_mode_cmd = messenger->DeclareMethod("evaluatedCrossSection", &H11BReaction::SetEvaluatedCrossSectionModeCommand, "Set evaluated total cross section");
  xs_mode_cmd.SetCandidates("tentori2023 model");

  messenger->DeclareMethod("crossSectionBiasFactor", &H11BReaction::SetCrossSectionBiasFactorCommand, "Set p + 11B cross-section sampling bias factor (>= 1)");

  messenger->DeclareMethod("enable675PrimaryAngularDistribution", &H11BReaction::SetEnable675PrimaryAngularDistributionCommand, "Enable 675-keV primary alpha A1/A2 angular distribution");
  messenger->DeclareMethod("675PrimaryAngularA1", &H11BReaction::Set675PrimaryAngularA1Command, "Set 675-keV primary alpha A1 coefficient");
  messenger->DeclareMethod("675PrimaryAngularA2", &H11BReaction::Set675PrimaryAngularA2Command, "Set 675-keV primary alpha A2 coefficient");

  messenger->DeclareMethod("enable675Alpha1SecondaryAngularCorrelation", &H11BReaction::SetEnable675Alpha1SecondaryAngularCorrelationCommand, "Enable 675-keV alpha1 secondary A2/A4 angular correlation");
  messenger->DeclareMethod("675Alpha1SecondaryA2", &H11BReaction::Set675Alpha1SecondaryA2Command, "Set 675-keV alpha1 secondary A2 coefficient");
  messenger->DeclareMethod("675Alpha1SecondaryA4", &H11BReaction::Set675Alpha1SecondaryA4Command, "Set 675-keV alpha1 secondary A4 coefficient");

  messenger->DeclareMethod("enable165PrimaryAngularDistribution", &H11BReaction::SetEnable165PrimaryAngularDistributionCommand, "Enable 165-keV primary alpha A1/A2 angular distribution");
  messenger->DeclareMethod("165PrimaryAngularA1", &H11BReaction::Set165PrimaryAngularA1Command, "Set 165-keV primary alpha A1 coefficient");
  messenger->DeclareMethod("165PrimaryAngularA2", &H11BReaction::Set165PrimaryAngularA2Command, "Set 165-keV primary alpha A2 coefficient");

  messenger->DeclareMethod("enable165Alpha1SecondaryAngularCorrelation", &H11BReaction::SetEnable165Alpha1SecondaryAngularCorrelationCommand, "Enable 165-keV alpha1 secondary A2/A4 angular correlation");
  messenger->DeclareMethod("165Alpha1SecondaryA2", &H11BReaction::Set165Alpha1SecondaryA2Command, "Set 165-keV alpha1 secondary A2 coefficient");
  messenger->DeclareMethod("165Alpha1SecondaryA4", &H11BReaction::Set165Alpha1SecondaryA4Command, "Set 165-keV alpha1 secondary A4 coefficient");

  messenger->DeclareMethod("enable165Gamma0AngularDistribution", &H11BReaction::SetEnable165Gamma0AngularDistributionCommand, "Enable 165-keV gamma0 A1/A2 angular distribution");
  messenger->DeclareMethod("165Gamma0AngularA1", &H11BReaction::Set165Gamma0AngularA1Command, "Set 165-keV gamma0 A1 coefficient");
  messenger->DeclareMethod("165Gamma0AngularA2", &H11BReaction::Set165Gamma0AngularA2Command, "Set 165-keV gamma0 A2 coefficient");

  messenger->DeclareMethod("enable165GammaCapture", &H11BReaction::Set165GammaCaptureEnabledCommand, "Enable or disable 165-keV gamma capture");
  messenger->DeclareMethod("enable675GammaCapture", &H11BReaction::Set675GammaCaptureEnabledCommand, "Enable or disable phenomenological 675-keV gamma capture");
  messenger->DeclareMethod("enable675GammaAngularDistribution", &H11BReaction::SetEnable675GammaAngularDistributionCommand, "Enable 675-keV gamma A1/A2 angular distribution");
  messenger->DeclareMethod("675GammaAngularA1", &H11BReaction::Set675GammaAngularA1Command, "Set 675-keV gamma A1 coefficient");
  messenger->DeclareMethod("675GammaAngularA2", &H11BReaction::Set675GammaAngularA2Command, "Set 675-keV gamma A2 coefficient");

  messenger->DeclareMethod("gammaBiasFactor", &H11BReaction::SetGammaBiasFactorCommand, "Set gamma-capture sampling bias factor (1 disables biasing)");
  messenger->DeclareMethod("printConfig", &H11BReaction::PrintConfigCommand, "Print current p + 11B runtime configuration");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4HadFinalState* H11BReaction::ApplyYourself(const G4HadProjectile& projectile, G4Nucleus& target)
{
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(G4HadFinalStateStatus::isAlive);

  G4int z_target = target.GetZ_asInt();
  const G4int a_target = target.GetA_asInt();
  const G4int a_project = projectile.GetDefinition()->GetBaryonNumber();

  if (a_project == 1 && z_target == 5 && a_target == 11 && SelectReactionChannel(projectile.GetKineticEnergy())) {
    theParticleChange.SetStatusChange(G4HadFinalStateStatus::stopAndKill);

    G4ParticleDefinition* target = G4IonTable::GetIonTable()->GetIon(z_target, a_target, 0. * CLHEP::keV);
    G4ParticleDefinition* product1 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition* product2 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition* product3 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);

    if (reaction_channel == ReactionChannel::GammaCapture12C) {
      GenerateGammaCapture(projectile, target);
    } else {
      ReactionKinematic(projectile, target, product1, product2, product3);
    }
  }

  return &theParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool H11BReaction::SelectReactionChannel(G4double kinetic_energy_lab)
{
  selected_components = H11BCrossSection::CalculateComponents(kinetic_energy_lab);
  if (selected_components.sigma_total_sampling_all <= 0.) return false;

  gamma_resonance = H11BGammaResonance::None;
  gamma_branch = H11BGammaBranch::None;

  const G4double sampled_cross_section = G4UniformRand() * selected_components.sigma_total_sampling_all;
  if (sampled_cross_section < selected_components.sigma_3alpha_sampling_total) {
    const G4double sampled_3alpha_cross_section = G4UniformRand() * selected_components.sigma_3alpha_sampling_total;
    if (sampled_3alpha_cross_section < selected_components.sigma_165_sampling) {
      reaction_channel = ReactionChannel::Resonance165;
      resonance_type = Resonance165;
    } else if (sampled_3alpha_cross_section < selected_components.sigma_165_sampling + selected_components.sigma_675_sampling) {
      reaction_channel = ReactionChannel::Resonance675;
      resonance_type = Resonance675;
    } else {
      reaction_channel = ReactionChannel::Background3Alpha;
      resonance_type = Resonance675;
    }
    return true;
  }

  const G4double sampled_gamma_cross_section = sampled_cross_section - selected_components.sigma_3alpha_sampling_total;
  reaction_channel = ReactionChannel::GammaCapture12C;

  if (sampled_gamma_cross_section < selected_components.sigma_gamma_165_0_sampling) {
    resonance_type = Resonance165;
    gamma_resonance = H11BGammaResonance::Resonance165;
    gamma_branch = H11BGammaBranch::Gamma165ToGroundState;
  } else if (sampled_gamma_cross_section < selected_components.sigma_gamma_165_0_sampling + selected_components.sigma_gamma_165_1_sampling) {
    resonance_type = Resonance165;
    gamma_resonance = H11BGammaResonance::Resonance165;
    gamma_branch = H11BGammaBranch::Gamma165ToFirstExcitedState;
  } else {
    resonance_type = Resonance675;
    gamma_resonance = H11BGammaResonance::Resonance675;
    gamma_branch = Select675GammaBranch();
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::ReactionKinematic(const G4HadProjectile& projectile, G4ParticleDefinition* target, G4ParticleDefinition* product1, G4ParticleDefinition* product2, G4ParticleDefinition* product3)
{
  reaction_data.Clear();

  G4double m_projectile = projectile.GetDefinition()->GetPDGMass();
  G4double m_target = target->GetPDGMass();

  G4ParticleDefinition* particle_4He = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
  G4double m_4He = particle_4He->GetPDGMass();

  G4ParticleDefinition* particle_8Be = G4IonTable::GetIonTable()->GetIon(4, 8, 0. * CLHEP::keV);
  G4double m_8Be = particle_8Be->GetPDGMass();

  G4double kinetic_lab_projectile = projectile.GetKineticEnergy(); // MeV
  G4double s = m_projectile * m_projectile + m_target * m_target + 2 * m_target * (m_projectile + kinetic_lab_projectile);
  G4double sqrt_s = std::sqrt(s);
  G4double ex_max = std::max(0., sqrt_s - m_4He - m_8Be);

  if (reaction_channel == ReactionChannel::Background3Alpha) {
    GenerateThreeBodyPhaseSpace(projectile, target, product1, product2, product3);
    return;
  }

  //
  // p+11B -> alpha+8Be
  G4double ex_8Be = 0.;
  G4bool is_alpha1_branch = false;
  G4bool is_alpha0_branch = false;
  if (ChooseAlpha0Channel())
    is_alpha0_branch = true;
  else
    is_alpha1_branch = true;

  if (is_alpha1_branch) {
    ex_8Be = Sample8Be2PlusExcitationEnergy(ex_max);
  }

  m_8Be += ex_8Be;

  G4double e_cm_projectile = (s + m_projectile * m_projectile - m_target * m_target) / (2 * sqrt_s);
  G4double e_cm_target = (s + m_target * m_target - m_projectile * m_projectile) / (2 * sqrt_s);
  G4double e_cm_alpha1 = (s + m_4He * m_4He - m_8Be * m_8Be) / (2 * sqrt_s);
  G4double e_cm_8Be = (s + m_8Be * m_8Be - m_4He * m_4He) / (2 * sqrt_s);
  G4double p_cm_projectile = SafeSqrt(e_cm_projectile * e_cm_projectile - m_projectile * m_projectile);
  G4double p_cm_target = SafeSqrt(e_cm_target * e_cm_target - m_target * m_target);
  G4double p_cm_alpha1 = SafeSqrt(e_cm_alpha1 * e_cm_alpha1 - m_4He * m_4He);
  G4double p_cm_8Be = SafeSqrt(e_cm_8Be * e_cm_8Be - m_8Be * m_8Be);

  G4LorentzVector lv_cm_total(0, 0, 0, sqrt_s);

  G4ThreeVector momentum_lab_projectile = projectile.GetMomentumDirection();
  momentum_lab_projectile *= SafeSqrt(kinetic_lab_projectile * kinetic_lab_projectile + 2 * kinetic_lab_projectile * m_projectile);
  G4ThreeVector momentum_lab_target(0, 0, 0);

  G4LorentzVector lv_lab_projectile(momentum_lab_projectile, m_projectile + kinetic_lab_projectile);
  G4LorentzVector lv_lab_target(momentum_lab_target, m_target);
  G4LorentzVector lv_lab_total = lv_lab_projectile + lv_lab_target;

  //
  // First breakup: 12C* -> primary alpha + 8Be.
  //
  // For the 165-keV and 675-keV resonances, optionally sample the
  // primary alpha direction relative to the beam axis with a fixed A1/A2
  // Legendre distribution.
  const G4double ecm_p11b_keV = (sqrt_s - m_projectile - m_target) / keV;
  G4int primary_angular_mode = 0; // 0: isotropic, 1: fixed A1/A2
  G4double primary_a1_used = 0.0;
  G4double primary_a2_used = 0.0;
  G4ThreeVector dir_cm_alpha1;
  if (resonance_type == Resonance165 && H11BConfig::GetEnable165PrimaryAngularDistribution()) {
    const auto primary_coefficients = H11BAngularDistribution::Get165PrimaryA1A2(ecm_p11b_keV, is_alpha1_branch);
    primary_a1_used = primary_coefficients.first;
    primary_a2_used = primary_coefficients.second;
    primary_angular_mode = 1;
    dir_cm_alpha1 = H11BAngularDistribution::Sample165PrimaryAlphaDirection(projectile.GetMomentumDirection(), is_alpha1_branch, ecm_p11b_keV);
  } else if (resonance_type == Resonance675 && H11BConfig::GetEnable675PrimaryAngularDistribution()) {
    primary_a1_used = H11BConfig::Get675PrimaryAngularA1();
    primary_a2_used = H11BConfig::Get675PrimaryAngularA2();
    primary_angular_mode = 1;
    dir_cm_alpha1 = H11BAngularDistribution::SampleDirectionFromLegendreA1A2(projectile.GetMomentumDirection(),
                                                                              primary_a1_used,
                                                                              primary_a2_used);
  } else {
    dir_cm_alpha1 = H11BAngularDistribution::SampleIsotropicDirection();
  }

  G4double theta_cm_alpha1 = dir_cm_alpha1.theta();
  G4double phi_cm_alpha1 = dir_cm_alpha1.phi();

  G4LorentzVector lv_cm_alpha1(p_cm_alpha1 * dir_cm_alpha1.x(), p_cm_alpha1 * dir_cm_alpha1.y(), p_cm_alpha1 * dir_cm_alpha1.z(), e_cm_alpha1);
  G4LorentzVector lv_cm_8Be = lv_cm_total - lv_cm_alpha1;

  G4LorentzVector lv_lab_alpha1 = lv_cm_alpha1;
  lv_lab_alpha1.boost(lv_lab_total.boostVector());

  G4LorentzVector lv_lab_8Be = lv_cm_8Be;
  lv_lab_8Be.boost(lv_lab_total.boostVector());
  const G4ThreeVector axis_8be_recoil_cm = lv_cm_8Be.vect().mag2() > 0.0 ? lv_cm_8Be.vect().unit() : -dir_cm_alpha1.unit();
#ifdef ReactionCout
  G4cout << "mass projectile " << m_projectile << G4endl;
  G4cout << "mass target " << m_target << G4endl;
  G4cout << "mass 4He " << m_4He << G4endl;
  G4cout << "mass 8Be " << m_8Be << G4endl;
  G4cout << "kinetic_lab_projectile : " << kinetic_lab_projectile << G4endl;
  G4cout << "ex_8Be : " << ex_8Be << " MeV" << G4endl;
  G4cout << "kinetic_energy_cm_projectile : " << e_cm_projectile - m_projectile << " MeV" << G4endl;
  G4cout << "kinetic_energy_cm_target : " << e_cm_target - m_target << " MeV" << G4endl;
  G4cout << "kinetic_energy_cm_alpha1 : " << e_cm_alpha1 - m_4He << " MeV" << G4endl;
  G4cout << "kinetic_energy_cm_8Be : " << e_cm_8Be - m_8Be << " MeV" << G4endl;
  G4cout << "total lab.px : " << lv_lab_total.px() << G4endl;
  G4cout << "total lab.py : " << lv_lab_total.py() << G4endl;
  G4cout << "total lab.pz : " << lv_lab_total.pz() << G4endl;
  G4cout << "total lab.kinetic : " << lv_lab_total.e() - m_projectile - m_target << G4endl;
  G4cout << "theta_cm_alpha1 " << theta_cm_alpha1 << G4endl;
  G4cout << "phi_cm_alpha1 " << phi_cm_alpha1 << G4endl;
  G4cout << "lv_cm_alpha1.px :        " << lv_cm_alpha1.px() << G4endl;
  G4cout << "lv_cm_alpha1.py :        " << lv_cm_alpha1.py() << G4endl;
  G4cout << "lv_cm_alpha1.pz :        " << lv_cm_alpha1.pz() << G4endl;
  G4cout << "lv_cm_alpha1.kinematic : " << lv_cm_alpha1.e() - m_4He << G4endl;
  G4cout << "lv_cm_8Be.px :        " << lv_cm_8Be.px() << G4endl;
  G4cout << "lv_cm_8Be.py :        " << lv_cm_8Be.py() << G4endl;
  G4cout << "lv_cm_8Be.pz :        " << lv_cm_8Be.pz() << G4endl;
  G4cout << "lv_cm_8Be.kinematic : " << lv_cm_8Be.e() - m_8Be << G4endl;
  G4cout << "lv_lab_alpha1.px :        " << lv_lab_alpha1.px() << G4endl;
  G4cout << "lv_lab_alpha1.py :        " << lv_lab_alpha1.py() << G4endl;
  G4cout << "lv_lab_alpha1.pz :        " << lv_lab_alpha1.pz() << G4endl;
  G4cout << "lv_lab_alpha1.kinematic : " << lv_lab_alpha1.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha1.theta :     " << lv_lab_alpha1.theta() << G4endl;
  G4cout << "lv_lab_alpha1.phi :       " << lv_lab_alpha1.phi() << G4endl;
  G4cout << "lv_lab_8Be.px :        " << lv_lab_8Be.px() << G4endl;
  G4cout << "lv_lab_8Be.py :        " << lv_lab_8Be.py() << G4endl;
  G4cout << "lv_lab_8Be.pz :        " << lv_lab_8Be.pz() << G4endl;
  G4cout << "lv_lab_8Be.kinematic : " << lv_lab_8Be.e() - m_8Be << G4endl;
  G4cout << "lv_lab_8Be.theta :     " << lv_lab_8Be.theta() << G4endl;
  G4cout << "lv_lab_8Be.phi :       " << lv_lab_8Be.phi() << G4endl;
#endif

  //
  // 8Be -> alpha+alpha
  G4double s_new = m_8Be * m_8Be;

  G4double sqrt_s_new = std::sqrt(s_new);
  G4double e_cm_alpha2 = sqrt_s_new / 2.;
  G4double e_cm_alpha3 = sqrt_s_new / 2.;
  G4double p_cm_alpha2 = SafeSqrt(e_cm_alpha2 * e_cm_alpha2 - m_4He * m_4He);
  G4double p_cm_alpha3 = SafeSqrt(e_cm_alpha3 * e_cm_alpha3 - m_4He * m_4He);

  G4ThreeVector dir_cm_alpha2;
  G4int h11b675_decay_model = 0;
  G4int secondary_angular_model = 0;
  G4double secondary_a2_used = 0.0;
  G4double secondary_a4_used = 0.0;
  G4double cos_theta_secondary_correlation = std::numeric_limits<G4double>::quiet_NaN();

  if (is_alpha1_branch && resonance_type == Resonance165) {
    // 165-keV alpha1 branch:
    //   12C*(16.11, 2+) -> alpha + 8Be(2+),
    //   8Be(2+) -> alpha + alpha.
    // The configurable Legendre model uses chi relative to the primary alpha
    // direction in the 12C CM frame.  Treado's fitted symmetry-axis shift is
    // intentionally ignored here.
    if (H11BConfig::GetEnable165Alpha1SecondaryAngularCorrelation()) {
      secondary_angular_model = 1;
      secondary_a2_used = H11BConfig::Get165Alpha1SecondaryA2();
      secondary_a4_used = H11BConfig::Get165Alpha1SecondaryA4();
      dir_cm_alpha2 = H11BAngularDistribution::SampleDirectionFromLegendreA2A4(
        lv_cm_alpha1.vect(), secondary_a2_used, secondary_a4_used);
      cos_theta_secondary_correlation = SafeCosBetween(lv_cm_alpha1.vect(), dir_cm_alpha2);
    } else {
      dir_cm_alpha2 = H11BAngularDistribution::SampleIsotropicDirection();
      cos_theta_secondary_correlation = SafeCosBetween(lv_cm_alpha1.vect(), dir_cm_alpha2);
    }
  } else if (is_alpha1_branch && resonance_type == Resonance675) {
    // 675-keV alpha1 branch:
    //   12C*(16.57/16.62, 2-) -> alpha + 8Be(2+)
    //   8Be(2+) -> alpha + alpha.
    // The secondary alpha correlation is controlled by a simple enable flag.
    if (H11BConfig::GetEnable675Alpha1SecondaryAngularCorrelation()) {
      h11b675_decay_model = 1;
      secondary_angular_model = 1;
      secondary_a2_used = H11BConfig::Get675Alpha1SecondaryA2();
      secondary_a4_used = H11BConfig::Get675Alpha1SecondaryA4();
      dir_cm_alpha2 = H11BAngularDistribution::Sample675Alpha1LegendreA2A4Direction(axis_8be_recoil_cm);
      cos_theta_secondary_correlation = SafeCosBetween(axis_8be_recoil_cm, dir_cm_alpha2);
    } else {
      h11b675_decay_model = 0;
      dir_cm_alpha2 = H11BAngularDistribution::SampleIsotropicDirection();
      cos_theta_secondary_correlation = SafeCosBetween(axis_8be_recoil_cm, dir_cm_alpha2);
    }
  } else {
    // alpha0 branch through 8Be(g.s.) or fallback: keep the old isotropic decay.
    dir_cm_alpha2 = H11BAngularDistribution::SampleIsotropicDirection();
  }
  G4double theta_cm_alpha2 = dir_cm_alpha2.theta();
  G4double phi_cm_alpha2 = dir_cm_alpha2.phi();
  const G4double cos_chi_8be_recoil = SafeCosBetween(axis_8be_recoil_cm, dir_cm_alpha2);
  const G4double phi_chi_8be_recoil = AzimuthAroundAxis(dir_cm_alpha2, axis_8be_recoil_cm);

  G4LorentzVector lv_cm_total_new(0., 0., 0., sqrt_s_new);
  G4LorentzVector lv_cm_alpha2(p_cm_alpha2 * dir_cm_alpha2.x(), p_cm_alpha2 * dir_cm_alpha2.y(), p_cm_alpha2 * dir_cm_alpha2.z(), e_cm_alpha2);
  G4LorentzVector lv_cm_alpha3 = lv_cm_total_new - lv_cm_alpha2;

  G4LorentzVector lv_lab_alpha2 = lv_cm_alpha2;
  lv_lab_alpha2.boost(lv_lab_8Be.boostVector());
  G4LorentzVector lv_lab_alpha3 = lv_cm_alpha3;
  lv_lab_alpha3.boost(lv_lab_8Be.boostVector());
#ifdef ReactionCout
  G4cout << "e_cm_alpha2 : " << e_cm_alpha2 - m_4He << G4endl;
  G4cout << "e_cm_alpha3 : " << e_cm_alpha3 - m_4He << G4endl;
  G4cout << "p_cm_alpha2 : " << p_cm_alpha2 << G4endl;
  G4cout << "p_cm_alpha3 : " << p_cm_alpha3 << G4endl;
  G4cout << "theta_cm_alpha2 " << theta_cm_alpha2 << G4endl;
  G4cout << "phi_cm_alpha2 " << phi_cm_alpha2 << G4endl;
  G4cout << "lv_cm_alpha2.px : " << lv_cm_alpha2.px() << G4endl;
  G4cout << "lv_cm_alpha2.py : " << lv_cm_alpha2.py() << G4endl;
  G4cout << "lv_cm_alpha2.pz : " << lv_cm_alpha2.pz() << G4endl;
  G4cout << "lv_cm_alpha2.kinematic : " << lv_cm_alpha2.e() - m_4He << G4endl;
  G4cout << "lv_cm_alpha3.px : " << lv_cm_alpha3.px() << G4endl;
  G4cout << "lv_cm_alpha3.py : " << lv_cm_alpha3.py() << G4endl;
  G4cout << "lv_cm_alpha3.pz : " << lv_cm_alpha3.pz() << G4endl;
  G4cout << "lv_cm_alpha3.kinematic : " << lv_cm_alpha3.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha2.px : " << lv_lab_alpha2.px() << G4endl;
  G4cout << "lv_lab_alpha2.py : " << lv_lab_alpha2.py() << G4endl;
  G4cout << "lv_lab_alpha2.pz : " << lv_lab_alpha2.pz() << G4endl;
  G4cout << "lv_lab_alpha2.kinematic : " << lv_lab_alpha2.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha2.theta :     " << lv_lab_alpha2.theta() << G4endl;
  G4cout << "lv_lab_alpha2.phi :       " << lv_lab_alpha2.phi() << G4endl;
  G4cout << "lv_lab_alpha3.px : " << lv_lab_alpha3.px() << G4endl;
  G4cout << "lv_lab_alpha3.py : " << lv_lab_alpha3.py() << G4endl;
  G4cout << "lv_lab_alpha3.pz : " << lv_lab_alpha3.pz() << G4endl;
  G4cout << "lv_lab_alpha3.kinematic : " << lv_lab_alpha3.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha3.theta :     " << lv_lab_alpha3.theta() << G4endl;
  G4cout << "lv_lab_alpha3.phi :       " << lv_lab_alpha3.phi() << G4endl;
#endif

  G4DynamicParticle* dynamic_product1 = new G4DynamicParticle(product1, lv_lab_alpha1);
  G4DynamicParticle* dynamic_product2 = new G4DynamicParticle(product2, lv_lab_alpha2);
  G4DynamicParticle* dynamic_product3 = new G4DynamicParticle(product3, lv_lab_alpha3);

  // Add the secondaries to the particle change stack
  // Changes being included in process
  theParticleChange.AddSecondary(dynamic_product1);
  theParticleChange.AddSecondary(dynamic_product2);
  theParticleChange.AddSecondary(dynamic_product3);

  // Build the three-alpha center-of-mass frame from the final-state alpha particles.
  // This is the proper frame for the Dalitz plot.
  G4LorentzVector lv_lab_3alpha = lv_lab_alpha1 + lv_lab_alpha2 + lv_lab_alpha3;

  G4ThreeVector beta_3alpha = lv_lab_3alpha.boostVector();

  G4LorentzVector lv_3alpha_cm_alpha1 = lv_lab_alpha1;
  G4LorentzVector lv_3alpha_cm_alpha2 = lv_lab_alpha2;
  G4LorentzVector lv_3alpha_cm_alpha3 = lv_lab_alpha3;

  lv_3alpha_cm_alpha1.boost(-beta_3alpha);
  lv_3alpha_cm_alpha2.boost(-beta_3alpha);
  lv_3alpha_cm_alpha3.boost(-beta_3alpha);

  // rootfile
  reaction_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  reaction_data.e_alpha1 = lv_lab_alpha1.e() - m_4He;
  reaction_data.e_alpha2 = lv_lab_alpha2.e() - m_4He;
  reaction_data.e_alpha3 = lv_lab_alpha3.e() - m_4He;
  // three-alpha CM-frame kinetic energies
  reaction_data.e_3alpha_cm_alpha1 = lv_3alpha_cm_alpha1.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha2 = lv_3alpha_cm_alpha2.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha3 = lv_3alpha_cm_alpha3.e() - m_4He;
  reaction_data.theta_lab_alpha1 = lv_lab_alpha1.theta();
  reaction_data.theta_lab_alpha2 = lv_lab_alpha2.theta();
  reaction_data.theta_lab_alpha3 = lv_lab_alpha3.theta();
  reaction_data.phi_lab_alpha1 = lv_lab_alpha1.phi();
  reaction_data.phi_lab_alpha2 = lv_lab_alpha2.phi();
  reaction_data.phi_lab_alpha3 = lv_lab_alpha3.phi();
  reaction_data.resonance_id = reaction_channel == ReactionChannel::Background3Alpha ? 0 : (resonance_type == Resonance675 ? 675 : 165);
  reaction_data.branch_id = is_alpha1_branch ? 1 : 0;
  reaction_data.reaction_channel = static_cast<G4int>(reaction_channel);
  reaction_data.background_mode = static_cast<G4int>(g_background_mode);
  reaction_data.e_cm_p11B = sqrt_s - m_projectile - m_target;
  reaction_data.projectile_kinetic_lab = kinetic_lab_projectile;
  reaction_data.projectile_px_lab = lv_lab_projectile.px();
  reaction_data.projectile_py_lab = lv_lab_projectile.py();
  reaction_data.projectile_pz_lab = lv_lab_projectile.pz();
  reaction_data.projectile_p_lab = lv_lab_projectile.vect().mag();
  reaction_data.projectile_theta_lab = lv_lab_projectile.vect().theta();
  reaction_data.projectile_phi_lab = lv_lab_projectile.vect().phi();
  reaction_data.ex_max_8Be = ex_max;
  reaction_data.eaa_8Be = ex_8Be + Ex8BeGroundAbove2Alpha;
  reaction_data.e_alpha8Be = std::max(0., ex_max - ex_8Be);
  const G4ThreeVector beam_axis = projectile.GetMomentumDirection();
  reaction_data.cos_theta_primary_cm = beam_axis.mag2() > 0. ? dir_cm_alpha1.unit().dot(beam_axis.unit()) : 0.;
  reaction_data.cos_chi_exit = SafeCosBetween(dir_cm_alpha1, dir_cm_alpha2);
  reaction_data.phi_primary_cm = phi_cm_alpha1;
  reaction_data.cos_chi_secondary_8be = cos_chi_8be_recoil;
  reaction_data.primary_angular_mode = primary_angular_mode;
  reaction_data.primary_a1_used = primary_a1_used;
  reaction_data.primary_a2_used = primary_a2_used;
  reaction_data.h11b675_decay_model = h11b675_decay_model;
  reaction_data.h11b675_decay_model_used = h11b675_decay_model;
  reaction_data.secondary_angular_model = secondary_angular_model;
  reaction_data.secondary_a2_used = secondary_a2_used;
  reaction_data.secondary_a4_used = secondary_a4_used;
  reaction_data.cos_theta_secondary_correlation = cos_theta_secondary_correlation;
  reaction_data.background_sequential_model = 0;
  reaction_data.cos_chi_675_internal = is_alpha1_branch && resonance_type == Resonance675 ? cos_chi_8be_recoil : 0.0;
  reaction_data.phi_chi_675_internal = is_alpha1_branch && resonance_type == Resonance675 ? phi_chi_8be_recoil : 0.0;
  reaction_data.opening_angle_alpha12_cm = SafeOpeningAngle(lv_3alpha_cm_alpha1.vect(), lv_3alpha_cm_alpha2.vect());
  reaction_data.opening_angle_alpha13_cm = SafeOpeningAngle(lv_3alpha_cm_alpha1.vect(), lv_3alpha_cm_alpha3.vect());
  reaction_data.opening_angle_alpha23_cm = SafeOpeningAngle(lv_3alpha_cm_alpha2.vect(), lv_3alpha_cm_alpha3.vect());
  reaction_data.e_alpha1_cm = reaction_data.e_3alpha_cm_alpha1;
  reaction_data.e_alpha2_cm = reaction_data.e_3alpha_cm_alpha2;
  reaction_data.e_alpha3_cm = reaction_data.e_3alpha_cm_alpha3;
  reaction_data.e_8be_excitation = ex_8Be;
  reaction_data.event_weight = ReactionChannelEventWeight(reaction_channel, selected_components);
  reaction_data.event_sampling_weight = reaction_data.event_weight;
  reaction_data.ex_8Be = ex_8Be;
  reaction_data.e_8Be = lv_lab_8Be.e() - m_8Be;
  reaction_data.theta_lab_8Be = lv_lab_8Be.theta();
  reaction_data.phi_lab_8Be = lv_lab_8Be.phi();

  G4ThreeVector reaction_position;
  auto event_manager = G4EventManager::GetEventManager();
  auto tracking_manager = event_manager ? event_manager->GetTrackingManager() : nullptr;
  auto current_track = tracking_manager ? tracking_manager->GetTrack() : nullptr;
  if (current_track) {
    reaction_position = current_track->GetPosition();
  }

  reaction_data.x = reaction_position.x();
  reaction_data.y = reaction_position.y();
  reaction_data.z = reaction_position.z();
  FillCrossSectionDiagnostics(reaction_data, selected_components);
  FillRuntimeConfigDiagnostics(reaction_data);
  reaction_data.sigma_eval_b = selected_components.sigma_eval / barn;
  reaction_data.sigma_165_model_b = selected_components.sigma_165_model / barn;
  reaction_data.sigma_675_model_b = selected_components.sigma_675_model / barn;
  reaction_data.sigma_165_used_b = selected_components.sigma_165 / barn;
  reaction_data.sigma_675_used_b = selected_components.sigma_675 / barn;
  reaction_data.sigma_total_used_b = selected_components.sigma_total / barn;
  reaction_data.sigma_model_sum_b = (selected_components.sigma_165_model + selected_components.sigma_675_model) / barn;
  reaction_data.model_scale_factor = reaction_data.sigma_model_sum_b > 0.0 ? (selected_components.sigma_165 + selected_components.sigma_675) / (selected_components.sigma_165_model + selected_components.sigma_675_model) : 1.0;
  reaction_data.sigma_background_b = selected_components.sigma_background / barn;
  if (selected_components.sigma_total > 0.) {
    reaction_data.channel_probability_165 = selected_components.sigma_165 / selected_components.sigma_total;
    reaction_data.channel_probability_675 = selected_components.sigma_675 / selected_components.sigma_total;
    reaction_data.channel_probability_background = selected_components.sigma_background / selected_components.sigma_total;
  }

  if (resonance_type == Resonance675 && is_alpha1_branch) {
    std::snprintf(reaction_data.reaction, sizeof(reaction_data.reaction), "%s_alpha1_%s", GetResonanceName(), H11BConfig::Get675Alpha1SecondaryAngularCorrelationName());
  } else if (is_alpha1_branch) {
    std::snprintf(reaction_data.reaction, sizeof(reaction_data.reaction), "%s_alpha1_%s", GetResonanceName(), H11BConfig::Get165Alpha1SecondaryAngularCorrelationName());
  } else if (is_alpha0_branch) {
    std::snprintf(reaction_data.reaction, sizeof(reaction_data.reaction), "%s_alpha0_isotropic", GetResonanceName());
  } else {
    // gamma
  }

  RunAction* run_action = static_cast<RunAction*>(const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
  run_action->GetRootIO()->FillReactionTree(reaction_data);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::GenerateGammaCapture(const G4HadProjectile& projectile, G4ParticleDefinition* target)
{
  reaction_data.Clear();

  const G4double m_projectile = projectile.GetDefinition()->GetPDGMass();
  const G4double m_target = target->GetPDGMass();
  const G4double kinetic_lab_projectile = projectile.GetKineticEnergy();
  const G4double s = m_projectile * m_projectile + m_target * m_target + 2 * m_target * (m_projectile + kinetic_lab_projectile);
  const G4double sqrt_s = std::sqrt(s);

  const G4ThreeVector momentum_lab_projectile = projectile.GetMomentumDirection() * SafeSqrt(kinetic_lab_projectile * kinetic_lab_projectile + 2 * kinetic_lab_projectile * m_projectile);
  const G4LorentzVector lv_lab_projectile(momentum_lab_projectile, m_projectile + kinetic_lab_projectile);
  const G4LorentzVector lv_lab_target(G4ThreeVector(0., 0., 0.), m_target);
  const G4LorentzVector lv_lab_total = lv_lab_projectile + lv_lab_target;
  const G4LorentzVector lv_cm_total(0., 0., 0., sqrt_s);

  G4ParticleDefinition* c12_ground = G4IonTable::GetIonTable()->GetIon(6, 12, 0. * CLHEP::keV);
  G4ParticleDefinition* gamma_particle = G4Gamma::Definition();

  reaction_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  reaction_data.resonance_id = static_cast<G4int>(gamma_resonance);
  reaction_data.branch_id = -1;
  reaction_data.reaction_channel = static_cast<G4int>(reaction_channel);
  reaction_data.background_mode = static_cast<G4int>(g_background_mode);
  reaction_data.gamma_resonance = static_cast<G4int>(gamma_resonance);
  reaction_data.gamma_branch = static_cast<G4int>(gamma_branch);
  reaction_data.e_cm_p11B = sqrt_s - m_projectile - m_target;
  reaction_data.projectile_kinetic_lab = kinetic_lab_projectile;
  reaction_data.projectile_px_lab = lv_lab_projectile.px();
  reaction_data.projectile_py_lab = lv_lab_projectile.py();
  reaction_data.projectile_pz_lab = lv_lab_projectile.pz();
  reaction_data.projectile_p_lab = lv_lab_projectile.vect().mag();
  reaction_data.projectile_theta_lab = lv_lab_projectile.vect().theta();
  reaction_data.projectile_phi_lab = lv_lab_projectile.vect().phi();
  reaction_data.gamma_bias_factor = H11BConfig::GetGammaBiasFactor();
  reaction_data.gamma_event_weight = 1.0;

  const G4double physical = GammaBranchPhysicalWeight(gamma_branch, selected_components);
  const G4double sampling = GammaBranchWeight(gamma_branch, selected_components);
  reaction_data.gamma_event_weight = sampling > 0.0 ? physical / sampling : 1.0;
  reaction_data.event_weight = reaction_data.gamma_event_weight;
  reaction_data.event_sampling_weight = reaction_data.gamma_event_weight;

  if (gamma_branch == H11BGammaBranch::Gamma165ToGroundState) {
    Generate165Gamma0(lv_cm_total, lv_lab_total, c12_ground, gamma_particle, projectile.GetMomentumDirection());
    std::snprintf(reaction_data.reaction,
                  sizeof(reaction_data.reaction),
                  "p11B_165_gamma0_%s",
                  H11BConfig::Get165Gamma0AngularDistributionName());
  } else if (gamma_branch == H11BGammaBranch::Gamma165ToFirstExcitedState) {
    Generate165Gamma1Cascade(lv_cm_total, lv_lab_total, c12_ground, gamma_particle);
    std::snprintf(reaction_data.reaction, sizeof(reaction_data.reaction), "p11B_165_gamma1_isotropicCascade");
  } else {
    Generate675GammaLine(lv_cm_total, lv_lab_total, c12_ground, gamma_particle, projectile.GetMomentumDirection());
    const H11BGammaLine* line = Find675GammaLine(gamma_branch);
    std::snprintf(reaction_data.reaction,
                  sizeof(reaction_data.reaction),
                  "p11B_675_gamma_relativeLineTable_%s_%s_scale1e-5",
                  line ? line->label : "unknown",
                  H11BConfig::Get675GammaAngularDistributionName());
  }

  G4ThreeVector reaction_position;
  auto event_manager = G4EventManager::GetEventManager();
  auto tracking_manager = event_manager ? event_manager->GetTrackingManager() : nullptr;
  auto current_track = tracking_manager ? tracking_manager->GetTrack() : nullptr;
  if (current_track) {
    reaction_position = current_track->GetPosition();
  }

  reaction_data.x = reaction_position.x();
  reaction_data.y = reaction_position.y();
  reaction_data.z = reaction_position.z();
  FillCrossSectionDiagnostics(reaction_data, selected_components);
  FillRuntimeConfigDiagnostics(reaction_data);

  RunAction* run_action = static_cast<RunAction*>(const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
  run_action->GetRootIO()->FillReactionTree(reaction_data);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Generate165Gamma0(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle, const G4ThreeVector& beam_axis_cm)
{
  const G4double m_final = c12_ground->GetPDGMass();
  const G4double m_initial = initial_cm.e();
  const G4double e_gamma = (m_initial * m_initial - m_final * m_final) / (2.0 * m_initial);
  const G4ThreeVector dir_gamma_cm =
    H11BConfig::GetEnable165Gamma0AngularDistribution()
      ? H11BAngularDistribution::Sample165Gamma0FixedA1A2Direction(beam_axis_cm)
      : H11BAngularDistribution::SampleIsotropicDirection();

  G4LorentzVector lv_cm_gamma(e_gamma * dir_gamma_cm, e_gamma);
  G4LorentzVector lv_cm_c12(-e_gamma * dir_gamma_cm, std::sqrt(m_final * m_final + e_gamma * e_gamma));

  G4LorentzVector lv_lab_gamma = lv_cm_gamma;
  G4LorentzVector lv_lab_c12 = lv_cm_c12;
  lv_lab_gamma.boost(initial_lab.boostVector());
  lv_lab_c12.boost(initial_lab.boostVector());

  theParticleChange.AddSecondary(new G4DynamicParticle(gamma_particle, lv_lab_gamma));
  theParticleChange.AddSecondary(new G4DynamicParticle(c12_ground, lv_lab_c12));

  reaction_data.gamma_angular_mode = static_cast<G4int>(H11BConfig::Get165Gamma0AngularModeForOutput());
  reaction_data.n_prompt_gammas = 1;
  reaction_data.gamma1_energy = lv_lab_gamma.e();
  reaction_data.gamma_primary_energy_MeV = e_gamma / MeV;
  reaction_data.gamma_final_state_energy_MeV = C12LevelGround / MeV;
  reaction_data.gamma_cascade_generated = 0;
  reaction_data.gamma1_theta_lab = lv_lab_gamma.theta();
  reaction_data.gamma1_phi_lab = lv_lab_gamma.phi();
  reaction_data.gamma1_theta_cm = lv_cm_gamma.theta();
  reaction_data.cos_theta_gamma_cm = SafeCosBetween(dir_gamma_cm, beam_axis_cm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Generate165Gamma1Cascade(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle)
{
  const G4double m_ground = c12_ground->GetPDGMass();
  const G4double m_intermediate = m_ground + C12Level4439;
  const G4double m_initial = initial_cm.e();

  const G4double e_gamma1 = (m_initial * m_initial - m_intermediate * m_intermediate) / (2.0 * m_initial);
  const G4ThreeVector dir_gamma1_cm = H11BAngularDistribution::SampleIsotropicDirection();
  G4LorentzVector lv_cm_gamma1(e_gamma1 * dir_gamma1_cm, e_gamma1);
  G4LorentzVector lv_cm_intermediate(-e_gamma1 * dir_gamma1_cm, std::sqrt(m_intermediate * m_intermediate + e_gamma1 * e_gamma1));

  G4LorentzVector lv_lab_gamma1 = lv_cm_gamma1;
  G4LorentzVector lv_lab_intermediate = lv_cm_intermediate;
  lv_lab_gamma1.boost(initial_lab.boostVector());
  lv_lab_intermediate.boost(initial_lab.boostVector());

  const G4double e_gamma2 = (m_intermediate * m_intermediate - m_ground * m_ground) / (2.0 * m_intermediate);
  const G4ThreeVector dir_gamma2_intermediate = H11BAngularDistribution::SampleIsotropicDirection();
  G4LorentzVector lv_intermediate_gamma2(e_gamma2 * dir_gamma2_intermediate, e_gamma2);
  G4LorentzVector lv_intermediate_c12(-e_gamma2 * dir_gamma2_intermediate, std::sqrt(m_ground * m_ground + e_gamma2 * e_gamma2));

  G4LorentzVector lv_lab_gamma2 = lv_intermediate_gamma2;
  G4LorentzVector lv_lab_c12 = lv_intermediate_c12;
  lv_lab_gamma2.boost(lv_lab_intermediate.boostVector());
  lv_lab_c12.boost(lv_lab_intermediate.boostVector());

  G4LorentzVector lv_cm_gamma2 = lv_intermediate_gamma2;
  lv_cm_gamma2.boost(lv_cm_intermediate.boostVector());

  theParticleChange.AddSecondary(new G4DynamicParticle(gamma_particle, lv_lab_gamma1));
  theParticleChange.AddSecondary(new G4DynamicParticle(gamma_particle, lv_lab_gamma2));
  theParticleChange.AddSecondary(new G4DynamicParticle(c12_ground, lv_lab_c12));

  reaction_data.gamma_angular_mode = static_cast<G4int>(H11BGammaAngularMode::Isotropic);
  reaction_data.n_prompt_gammas = 2;
  reaction_data.gamma1_energy = lv_lab_gamma1.e();
  reaction_data.gamma2_energy = lv_lab_gamma2.e();
  reaction_data.gamma_primary_energy_MeV = e_gamma1 / MeV;
  reaction_data.gamma_final_state_energy_MeV = C12Level4439 / MeV;
  reaction_data.gamma_cascade_generated = 1;
  reaction_data.gamma1_theta_lab = lv_lab_gamma1.theta();
  reaction_data.gamma2_theta_lab = lv_lab_gamma2.theta();
  reaction_data.gamma1_phi_lab = lv_lab_gamma1.phi();
  reaction_data.gamma2_phi_lab = lv_lab_gamma2.phi();
  reaction_data.gamma1_theta_cm = lv_cm_gamma1.theta();
  reaction_data.gamma2_theta_cm = lv_cm_gamma2.theta();
  reaction_data.cos_theta_gamma_cm = std::numeric_limits<G4double>::quiet_NaN();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::Generate675GammaLine(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle, const G4ThreeVector& beam_axis_cm)
{
  const H11BGammaLine* line = Find675GammaLine(gamma_branch);
  const G4double final_state_energy = line ? line->final_state_energy : C12Level4439;
  G4ParticleDefinition* c12_final = final_state_energy > 0.0 ? G4IonTable::GetIonTable()->GetIon(6, 12, final_state_energy) : c12_ground;

  const G4double m_final = c12_final->GetPDGMass();
  const G4double m_initial = initial_cm.e();
  const G4double e_gamma = (m_initial * m_initial - m_final * m_final) / (2.0 * m_initial);
  const G4ThreeVector dir_gamma_cm =
    H11BConfig::GetEnable675GammaAngularDistribution()
      ? H11BAngularDistribution::SampleDirectionFromLegendreA1A2(beam_axis_cm,
                                                                  H11BConfig::Get675GammaAngularA1(),
                                                                  H11BConfig::Get675GammaAngularA2())
      : H11BAngularDistribution::SampleIsotropicDirection();

  G4LorentzVector lv_cm_gamma(e_gamma * dir_gamma_cm, e_gamma);
  G4LorentzVector lv_cm_c12_final(-e_gamma * dir_gamma_cm, std::sqrt(m_final * m_final + e_gamma * e_gamma));

  G4LorentzVector lv_lab_gamma = lv_cm_gamma;
  G4LorentzVector lv_lab_c12_final = lv_cm_c12_final;
  lv_lab_gamma.boost(initial_lab.boostVector());
  lv_lab_c12_final.boost(initial_lab.boostVector());

  theParticleChange.AddSecondary(new G4DynamicParticle(gamma_particle, lv_lab_gamma));

  reaction_data.gamma_angular_mode = static_cast<G4int>(H11BConfig::Get675GammaAngularModeForOutput());
  reaction_data.n_prompt_gammas = 1;
  reaction_data.gamma1_energy = lv_lab_gamma.e();
  reaction_data.gamma_primary_energy_MeV = e_gamma / MeV;
  reaction_data.gamma_final_state_energy_MeV = final_state_energy / MeV;
  reaction_data.gamma_relative_intensity_used = line ? line->relative_intensity : 100.0;
  const G4double intensity_sum = Sum675GammaRelativeIntensity();
  reaction_data.gamma_branch_fraction_used =
    intensity_sum > 0.0 && line ? line->relative_intensity / intensity_sum : 1.0;
  reaction_data.gamma_branch_is_upper_limit = line && line->is_upper_limit ? 1 : 0;
  reaction_data.gamma_branch_from_relative_table = 1;
  reaction_data.gamma1_theta_lab = lv_lab_gamma.theta();
  reaction_data.gamma1_phi_lab = lv_lab_gamma.phi();
  reaction_data.gamma1_theta_cm = lv_cm_gamma.theta();
  reaction_data.cos_theta_gamma_cm = SafeCosBetween(dir_gamma_cm, beam_axis_cm);

  if (gamma_branch == H11BGammaBranch::Gamma675To4439State) {
    const G4double m_ground = c12_ground->GetPDGMass();
    const G4double m_intermediate = c12_final->GetPDGMass();
    const G4double e_gamma2 = (m_intermediate * m_intermediate - m_ground * m_ground) / (2.0 * m_intermediate);
    const G4ThreeVector dir_gamma2_intermediate = H11BAngularDistribution::SampleIsotropicDirection();
    G4LorentzVector lv_intermediate_gamma2(e_gamma2 * dir_gamma2_intermediate, e_gamma2);
    G4LorentzVector lv_intermediate_c12(-e_gamma2 * dir_gamma2_intermediate, std::sqrt(m_ground * m_ground + e_gamma2 * e_gamma2));

    G4LorentzVector lv_lab_gamma2 = lv_intermediate_gamma2;
    G4LorentzVector lv_lab_c12 = lv_intermediate_c12;
    lv_lab_gamma2.boost(lv_lab_c12_final.boostVector());
    lv_lab_c12.boost(lv_lab_c12_final.boostVector());

    G4LorentzVector lv_cm_gamma2 = lv_intermediate_gamma2;
    lv_cm_gamma2.boost(lv_cm_c12_final.boostVector());

    theParticleChange.AddSecondary(new G4DynamicParticle(gamma_particle, lv_lab_gamma2));
    theParticleChange.AddSecondary(new G4DynamicParticle(c12_ground, lv_lab_c12));

    reaction_data.n_prompt_gammas = 2;
    reaction_data.gamma2_energy = lv_lab_gamma2.e();
    reaction_data.gamma2_theta_lab = lv_lab_gamma2.theta();
    reaction_data.gamma2_phi_lab = lv_lab_gamma2.phi();
    reaction_data.gamma2_theta_cm = lv_cm_gamma2.theta();
    reaction_data.gamma_cascade_generated = 1;
  } else {
    theParticleChange.AddSecondary(new G4DynamicParticle(c12_final, lv_lab_c12_final));
    reaction_data.gamma_cascade_generated = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H11BReaction::GenerateThreeBodyPhaseSpace(const G4HadProjectile& projectile, G4ParticleDefinition* target, G4ParticleDefinition* product1, G4ParticleDefinition* product2, G4ParticleDefinition* product3)
{
  reaction_data.Clear();

  const G4double m_projectile = projectile.GetDefinition()->GetPDGMass();
  const G4double m_target = target->GetPDGMass();
  const G4double m_4He = product1->GetPDGMass();

  const G4double kinetic_lab_projectile = projectile.GetKineticEnergy();
  const G4double s = m_projectile * m_projectile + m_target * m_target + 2 * m_target * (m_projectile + kinetic_lab_projectile);
  const G4double sqrt_s = std::sqrt(s);
  const G4double ex_max = std::max(0., sqrt_s - 3.0 * m_4He);

  const G4double p_lab_projectile = SafeSqrt(kinetic_lab_projectile * kinetic_lab_projectile + 2 * kinetic_lab_projectile * m_projectile);
  G4ThreeVector momentum_lab_projectile = projectile.GetMomentumDirection() * p_lab_projectile;
  G4LorentzVector lv_lab_projectile(momentum_lab_projectile, m_projectile + kinetic_lab_projectile);
  G4LorentzVector lv_lab_target(G4ThreeVector(0., 0., 0.), m_target);
  G4LorentzVector lv_lab_total = lv_lab_projectile + lv_lab_target;

  const G4double m23 = SamplePhaseSpaceM23(sqrt_s, m_4He);
  const G4double e_cm_alpha1 = (s + m_4He * m_4He - m23 * m23) / (2.0 * sqrt_s);
  const G4double p_cm_alpha1 = SafeSqrt(e_cm_alpha1 * e_cm_alpha1 - m_4He * m_4He);

  const G4ThreeVector dir_cm_alpha1 = H11BAngularDistribution::SampleIsotropicDirection();
  G4LorentzVector lv_cm_total(0., 0., 0., sqrt_s);
  G4LorentzVector lv_cm_alpha1(p_cm_alpha1 * dir_cm_alpha1.x(), p_cm_alpha1 * dir_cm_alpha1.y(), p_cm_alpha1 * dir_cm_alpha1.z(), e_cm_alpha1);
  G4LorentzVector lv_cm_23 = lv_cm_total - lv_cm_alpha1;

  const G4double e_23_alpha2 = m23 / 2.0;
  const G4double p_23_alpha2 = SafeSqrt(e_23_alpha2 * e_23_alpha2 - m_4He * m_4He);
  const G4ThreeVector dir_23_alpha2 = H11BAngularDistribution::SampleIsotropicDirection();

  G4LorentzVector lv_23_alpha2(p_23_alpha2 * dir_23_alpha2.x(), p_23_alpha2 * dir_23_alpha2.y(), p_23_alpha2 * dir_23_alpha2.z(), e_23_alpha2);
  G4LorentzVector lv_23_alpha3(-lv_23_alpha2.vect(), e_23_alpha2);

  G4LorentzVector lv_cm_alpha2 = lv_23_alpha2;
  lv_cm_alpha2.boost(lv_cm_23.boostVector());
  G4LorentzVector lv_cm_alpha3 = lv_23_alpha3;
  lv_cm_alpha3.boost(lv_cm_23.boostVector());

  G4LorentzVector lv_lab_alpha1 = lv_cm_alpha1;
  lv_lab_alpha1.boost(lv_lab_total.boostVector());
  G4LorentzVector lv_lab_alpha2 = lv_cm_alpha2;
  lv_lab_alpha2.boost(lv_lab_total.boostVector());
  G4LorentzVector lv_lab_alpha3 = lv_cm_alpha3;
  lv_lab_alpha3.boost(lv_lab_total.boostVector());

  theParticleChange.AddSecondary(new G4DynamicParticle(product1, lv_lab_alpha1));
  theParticleChange.AddSecondary(new G4DynamicParticle(product2, lv_lab_alpha2));
  theParticleChange.AddSecondary(new G4DynamicParticle(product3, lv_lab_alpha3));

  G4LorentzVector lv_lab_3alpha = lv_lab_alpha1 + lv_lab_alpha2 + lv_lab_alpha3;
  const G4ThreeVector beta_3alpha = lv_lab_3alpha.boostVector();
  G4LorentzVector lv_3alpha_cm_alpha1 = lv_lab_alpha1;
  G4LorentzVector lv_3alpha_cm_alpha2 = lv_lab_alpha2;
  G4LorentzVector lv_3alpha_cm_alpha3 = lv_lab_alpha3;
  lv_3alpha_cm_alpha1.boost(-beta_3alpha);
  lv_3alpha_cm_alpha2.boost(-beta_3alpha);
  lv_3alpha_cm_alpha3.boost(-beta_3alpha);

  reaction_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  reaction_data.e_alpha1 = lv_lab_alpha1.e() - m_4He;
  reaction_data.e_alpha2 = lv_lab_alpha2.e() - m_4He;
  reaction_data.e_alpha3 = lv_lab_alpha3.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha1 = lv_3alpha_cm_alpha1.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha2 = lv_3alpha_cm_alpha2.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha3 = lv_3alpha_cm_alpha3.e() - m_4He;
  reaction_data.theta_lab_alpha1 = lv_lab_alpha1.theta();
  reaction_data.theta_lab_alpha2 = lv_lab_alpha2.theta();
  reaction_data.theta_lab_alpha3 = lv_lab_alpha3.theta();
  reaction_data.phi_lab_alpha1 = lv_lab_alpha1.phi();
  reaction_data.phi_lab_alpha2 = lv_lab_alpha2.phi();
  reaction_data.phi_lab_alpha3 = lv_lab_alpha3.phi();
  reaction_data.resonance_id = 0;
  reaction_data.branch_id = -1;
  reaction_data.reaction_channel = static_cast<G4int>(reaction_channel);
  reaction_data.background_mode = static_cast<G4int>(g_background_mode);
  reaction_data.e_cm_p11B = sqrt_s - m_projectile - m_target;
  reaction_data.ex_max_8Be = ex_max;
  reaction_data.eaa_8Be = 0.;
  reaction_data.e_alpha8Be = 0.;
  reaction_data.projectile_kinetic_lab = kinetic_lab_projectile;
  reaction_data.projectile_px_lab = lv_lab_projectile.px();
  reaction_data.projectile_py_lab = lv_lab_projectile.py();
  reaction_data.projectile_pz_lab = lv_lab_projectile.pz();
  reaction_data.projectile_p_lab = lv_lab_projectile.vect().mag();
  reaction_data.projectile_theta_lab = lv_lab_projectile.vect().theta();
  reaction_data.projectile_phi_lab = lv_lab_projectile.vect().phi();

  const G4ThreeVector beam_axis = projectile.GetMomentumDirection();
  reaction_data.cos_theta_primary_cm = beam_axis.mag2() > 0. && lv_cm_alpha1.vect().mag2() > 0. ? lv_cm_alpha1.vect().unit().dot(beam_axis.unit()) : 0.;
  reaction_data.cos_chi_exit = lv_cm_alpha1.vect().mag2() > 0. && lv_cm_alpha2.vect().mag2() > 0. ? lv_cm_alpha1.vect().unit().dot(lv_cm_alpha2.vect().unit()) : 0.;
  reaction_data.phi_primary_cm = lv_cm_alpha1.phi();
  reaction_data.cos_chi_secondary_8be = reaction_data.cos_chi_exit;
  reaction_data.primary_angular_mode = 0;
  reaction_data.primary_a1_used = 0.0;
  reaction_data.primary_a2_used = 0.0;
  reaction_data.h11b675_decay_model = 0;
  reaction_data.h11b675_decay_model_used = 0;
  reaction_data.background_sequential_model = 0;
  reaction_data.cos_chi_675_internal = 0.0;
  reaction_data.phi_chi_675_internal = 0.0;
  reaction_data.opening_angle_alpha12_cm = SafeOpeningAngle(lv_3alpha_cm_alpha1.vect(), lv_3alpha_cm_alpha2.vect());
  reaction_data.opening_angle_alpha13_cm = SafeOpeningAngle(lv_3alpha_cm_alpha1.vect(), lv_3alpha_cm_alpha3.vect());
  reaction_data.opening_angle_alpha23_cm = SafeOpeningAngle(lv_3alpha_cm_alpha2.vect(), lv_3alpha_cm_alpha3.vect());
  reaction_data.e_alpha1_cm = reaction_data.e_3alpha_cm_alpha1;
  reaction_data.e_alpha2_cm = reaction_data.e_3alpha_cm_alpha2;
  reaction_data.e_alpha3_cm = reaction_data.e_3alpha_cm_alpha3;
  reaction_data.e_8be_excitation = 0.0;
  reaction_data.event_weight = ReactionChannelEventWeight(reaction_channel, selected_components);
  reaction_data.event_sampling_weight = reaction_data.event_weight;

  G4ThreeVector reaction_position;
  auto event_manager = G4EventManager::GetEventManager();
  auto tracking_manager = event_manager ? event_manager->GetTrackingManager() : nullptr;
  auto current_track = tracking_manager ? tracking_manager->GetTrack() : nullptr;
  if (current_track) {
    reaction_position = current_track->GetPosition();
  }

  reaction_data.x = reaction_position.x();
  reaction_data.y = reaction_position.y();
  reaction_data.z = reaction_position.z();
  FillCrossSectionDiagnostics(reaction_data, selected_components);
  FillRuntimeConfigDiagnostics(reaction_data);
  reaction_data.sigma_eval_b = selected_components.sigma_eval / barn;
  reaction_data.sigma_165_model_b = selected_components.sigma_165_model / barn;
  reaction_data.sigma_675_model_b = selected_components.sigma_675_model / barn;
  reaction_data.sigma_165_used_b = selected_components.sigma_165 / barn;
  reaction_data.sigma_675_used_b = selected_components.sigma_675 / barn;
  reaction_data.sigma_total_used_b = selected_components.sigma_total / barn;
  reaction_data.sigma_model_sum_b = (selected_components.sigma_165_model + selected_components.sigma_675_model) / barn;
  reaction_data.model_scale_factor = reaction_data.sigma_model_sum_b > 0.0 ? (selected_components.sigma_165 + selected_components.sigma_675) / (selected_components.sigma_165_model + selected_components.sigma_675_model) : 1.0;
  reaction_data.sigma_background_b = selected_components.sigma_background / barn;
  if (selected_components.sigma_total > 0.) {
    reaction_data.channel_probability_165 = selected_components.sigma_165 / selected_components.sigma_total;
    reaction_data.channel_probability_675 = selected_components.sigma_675 / selected_components.sigma_total;
    reaction_data.channel_probability_background = selected_components.sigma_background / selected_components.sigma_total;
  }

  std::snprintf(reaction_data.reaction, sizeof(reaction_data.reaction), "background3alpha_phaseSpace");

  RunAction* run_action = static_cast<RunAction*>(const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
  run_action->GetRootIO()->FillReactionTree(reaction_data);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BReaction::SamplePhaseSpaceM23(G4double sqrt_s, G4double m_alpha) const
{
  const G4double min_m23 = 2.0 * m_alpha;
  const G4double max_m23 = sqrt_s - m_alpha;
  if (max_m23 <= min_m23) return min_m23;

  auto two_body_momentum = [](G4double parent_m, G4double m_a, G4double m_b) {
    if (parent_m <= m_a + m_b) return 0.0;

    const G4double a = parent_m * parent_m - (m_a + m_b) * (m_a + m_b);
    const G4double b = parent_m * parent_m - (m_a - m_b) * (m_a - m_b);
    return SafeSqrt(a * b) / (2.0 * parent_m);
  };

  auto weight = [&](G4double m23) {
    return two_body_momentum(sqrt_s, m_alpha, m23) * two_body_momentum(m23, m_alpha, m_alpha);
  };

  G4double max_weight = 0.;
  constexpr G4int n_scan = 128;
  for (G4int i = 0; i <= n_scan; ++i) {
    const G4double x = min_m23 + (max_m23 - min_m23) * static_cast<G4double>(i) / n_scan;
    max_weight = std::max(max_weight, weight(x));
  }

  if (max_weight <= 0.) return min_m23;

  for (G4int i = 0; i < 10000; ++i) {
    const G4double candidate = min_m23 + (max_m23 - min_m23) * G4UniformRand();
    if (G4UniformRand() * max_weight <= weight(candidate)) return candidate;
  }

  return 0.5 * (min_m23 + max_m23);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BReaction::Sample8Be2PlusExcitationEnergy(G4double ex_max)
{
  if (ex_max <= 0.) return 0.;

  // The CDF is cached for each kinematic upper limit and resonance type.
  // The resonance type matters because the first breakup
  //   12C* -> alpha + 8Be(2+)
  // uses different alpha+8Be penetrabilities: L=2 for 165 keV and
  // L=1/L=3 mixture for 675 keV.
  static thread_local std::map<G4int, ExcitationCdfTable> cdf_cache;
  G4int ex_key = static_cast<G4int>(10.0 * ex_max / keV + 0.5);
  G4int key = 100 * ex_key + 10 * (resonance_type == Resonance675 ? 1 : 0);

  auto table_it = cdf_cache.find(key);
  if (table_it == cdf_cache.end()) {
    const G4double step = 5. * keV;
    const G4int n_points = std::max(32, static_cast<G4int>(ex_max / step) + 1);

    ExcitationCdfTable table;
    table.ex_max = ex_max;
    table.x.resize(n_points);
    table.cdf.resize(n_points, 0.);

    std::vector<G4double> weights(n_points, 0.);
    for (G4int i = 0; i < n_points; ++i) {
      table.x[i] = static_cast<G4double>(i) * ex_max / (n_points - 1);
      weights[i] = std::max(0., Weight8Be2Plus(table.x[i], ex_max));
    }

    for (G4int i = 1; i < n_points; ++i) {
      const G4double dx = table.x[i] - table.x[i - 1];
      const G4double area = 0.5 * (weights[i - 1] + weights[i]) * dx;
      table.cdf[i] = table.cdf[i - 1] + area;
    }

    const G4double total_area = table.cdf.back();
    if (total_area > 0.) {
      for (auto& value : table.cdf) {
        value /= total_area;
      }
    }

    table_it = cdf_cache.emplace(key, std::move(table)).first;
  }

  const auto& table = table_it->second;
  if (table.x.empty() || table.cdf.empty() || table.cdf.back() <= 0.) return 0.;

  const G4double u = G4UniformRand();
  const auto upper = std::lower_bound(table.cdf.begin(), table.cdf.end(), u);
  if (upper == table.cdf.begin()) return 0.;
  if (upper == table.cdf.end()) return ex_max;

  const auto i = static_cast<G4int>(std::distance(table.cdf.begin(), upper));
  const G4double cdf0 = table.cdf[i - 1];
  const G4double cdf1 = table.cdf[i];
  const G4double fraction = cdf1 > cdf0 ? (u - cdf0) / (cdf1 - cdf0) : 0.;
  const G4double sampled_ex = table.x[i - 1] + fraction * (table.x[i] - table.x[i - 1]);

  return std::min(sampled_ex, ex_max);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BReaction::Weight8Be2Plus(G4double ex_8Be, G4double ex_max) const
{
  // ex_8Be is relative to the 8Be ground state, which is the quantity added
  // to the 8Be mass in ReactionKinematic(). For the alpha-alpha line shape,
  // the relevant energy is instead the relative energy above the 2-alpha
  // threshold:
  //     Eaa = ex_8Be + Ex8BeGroundAbove2Alpha.
  G4double eaa = ex_8Be + Ex8BeGroundAbove2Alpha;
  if (eaa <= 0.) return 0.;

  // Energy-dependent width for 8Be(2+) -> alpha + alpha.
  // Default: full Coulomb penetrability ratio
  //     Gamma(E) = Gamma_R * P_2(E,a) / P_2(E_R,a),
  // where P_L = rho/[F_L(eta,rho)^2 + G_L(eta,rho)^2].
  // If the full-Coulomb flag is disabled, the code falls back to the older
  // threshold approximation Gamma(E)=Gamma_R*(E/E_R)^(5/2).
  G4double gamma_e = Gamma8Be2Plus * std::pow(eaa / Eaa8Be2Plus, 2.5);
  if (H11BUseFullCoulombPenetrabilityFor8Be2Plus) {
    const G4double ratio = CoulombPenetrability::AlphaAlphaL2Ratio(eaa / keV, Eaa8Be2Plus / keV);
    gamma_e = Gamma8Be2Plus * ratio;
  }

  G4double delta_e = eaa - Eaa8Be2Plus;
  G4double w_8Be_line_shape = gamma_e / (delta_e * delta_e + gamma_e * gamma_e / 4.);

  // First-breakup penetrability for
  //   12C* -> alpha + 8Be(2+).
  // This is different from the alpha-alpha penetrability above.  Increasing
  // ex_8Be reduces the alpha+8Be relative energy and must therefore be
  // suppressed by the corresponding exit-channel penetrability.
  G4double e_alpha8Be = ex_max - ex_8Be;
  G4double e_alpha8Be_ref = ex_max - Ex8Be2Plus;
  G4double w_first_breakup = WeightAlpha8BeFirstBreakup(e_alpha8Be, e_alpha8Be_ref);

  return w_first_breakup * w_8Be_line_shape;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BReaction::WeightAlpha8BeFirstBreakup(G4double e_alpha8Be, G4double e_alpha8Be_ref) const
{
  if (e_alpha8Be <= 0.) return 0.;

  // If the nominal 8Be(2+) centroid is outside the available phase space,
  // do not renormalize to an unphysical reference. This should not occur for
  // the 165-keV or 675-keV resonances, but keeps the function robust.
  if (e_alpha8Be_ref <= 0.) return 1.;

  if (H11BUseFullCoulombPenetrabilityForAlpha8BeFirstBreakup) {
    if (resonance_type == Resonance165) {
      return CoulombPenetrability::Alpha8BeRatio(e_alpha8Be / keV, e_alpha8Be_ref / keV, H11B165Alpha1ExitOrbitalL);
    }

    // For the 675-keV 2- alpha1 branch, use the retained L=1/L=3
    // first-breakup penetrability mixture.  The secondary A2/A4 angular
    // correlation affects only the later 8Be(2+) -> alpha + alpha direction.
    const G4double f1 = H11B675ExitL1Fraction;
    return f1 *
               CoulombPenetrability::Alpha8BeRatio(e_alpha8Be / keV, e_alpha8Be_ref / keV, H11B675Alpha1ExitOrbitalL1) +
           (1.0 - f1) *
               CoulombPenetrability::Alpha8BeRatio(e_alpha8Be / keV, e_alpha8Be_ref / keV, H11B675Alpha1ExitOrbitalL3);
  }

  // Fallback without full Coulomb functions: threshold power law
  // P_L(E) ~ E^(L+1/2), equivalent to q^(2L+1).
  if (resonance_type == Resonance165) {
    return std::pow(e_alpha8Be / e_alpha8Be_ref, H11B165Alpha1ExitOrbitalL + 0.5);
  }

  const G4double f1 = H11B675ExitL1Fraction;
  return f1 * std::pow(e_alpha8Be / e_alpha8Be_ref, H11B675Alpha1ExitOrbitalL1 + 0.5) +
         (1.0 - f1) * std::pow(e_alpha8Be / e_alpha8Be_ref, H11B675Alpha1ExitOrbitalL3 + 0.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool H11BReaction::ChooseAlpha0Channel() const
{
  // The relative alpha0/alpha1 yield is taken directly from the partial widths
  // of the selected resonance. alpha0 means 8Be ground state; alpha1 means
  // 8Be(2+) excited state.
  G4double probability_alpha0 = GetAlpha0BranchingRatio();
  G4double rand = G4UniformRand();

  if (rand < probability_alpha0) return true; // alpha0 -> 8Be ground state
  return false;                               // alpha1 -> 8Be excited state
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double H11BReaction::GetAlpha0BranchingRatio() const
{
  G4double alpha0_width = H11B165Alpha0Width / keV;
  G4double alpha1_width = H11B165Alpha1Width / keV;

  if (resonance_type == Resonance675) {
    alpha0_width = H11B675Alpha0Width / keV;
    alpha1_width = H11B675Alpha1Width / keV;
  }

  G4double total_alpha_width = alpha0_width + alpha1_width;
  if (total_alpha_width <= 0.) return 0.;

  return alpha0_width / total_alpha_width;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const char* H11BReaction::GetResonanceName() const
{
  if (resonance_type == Resonance675) return "p11B_675";
  return "p11B_165";
}
