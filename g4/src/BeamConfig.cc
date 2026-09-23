#include "BeamConfig.hh"

#include "Constants.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cctype>

namespace {
BeamProfile g_profile = BeamProfile::UniformDisk;

// Preserve the previous normal-beam behavior by default.
G4double g_radius = BeamR;
G4double g_sigma_x = 1. * mm;
G4double g_sigma_y = 1. * mm;
G4double g_offset_x = 0.;
G4double g_offset_y = 0.;
G4double g_z = BeamZ;

G4String NormalizeProfileName(G4String profile)
{
  std::transform(profile.begin(),
                 profile.end(),
                 profile.begin(),
                 [](unsigned char value) {
                   return static_cast<char>(std::tolower(value));
                 });

  profile.erase(std::remove_if(profile.begin(),
                               profile.end(),
                               [](unsigned char value) {
                                 return std::isspace(value) != 0
                                     || value == '-'
                                     || value == '_';
                               }),
                profile.end());

  return profile;
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BeamConfig::BeamConfig()
{
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BeamConfig::~BeamConfig() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
BeamProfile BeamConfig::GetProfile()
{
  return g_profile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String BeamConfig::GetProfileName()
{
  switch (g_profile) {
  case BeamProfile::Point:
    return "point";
  case BeamProfile::UniformDisk:
    return "uniformDisk";
  case BeamProfile::Gaussian:
    return "gaussian";
  }

  return "unknown";
}

G4double BeamConfig::GetRadius() { return g_radius; }
G4double BeamConfig::GetSigmaX() { return g_sigma_x; }
G4double BeamConfig::GetSigmaY() { return g_sigma_y; }
G4double BeamConfig::GetOffsetX() { return g_offset_x; }
G4double BeamConfig::GetOffsetY() { return g_offset_y; }
G4double BeamConfig::GetZ() { return g_z; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::SetProfile(const G4String& profile)
{
  const G4String normalized = NormalizeProfileName(profile);

  if (normalized == "point") {
    g_profile = BeamProfile::Point;
    return;
  }

  if (normalized == "uniformdisk" || normalized == "disk") {
    g_profile = BeamProfile::UniformDisk;
    return;
  }

  if (normalized == "gaussian" || normalized == "gauss") {
    g_profile = BeamProfile::Gaussian;
    return;
  }

  G4cerr << "BeamConfig: unsupported profile \"" << profile << "\". "
         << "Valid values are point, uniformDisk and gaussian. "
         << "Keeping profile = " << GetProfileName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::SetRadius(G4double radius)
{
  if (radius < 0.) {
    G4cerr << "BeamConfig: radius must be >= 0. Keeping "
           << g_radius / mm << " mm." << G4endl;
    return;
  }

  g_radius = radius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::SetSigmaX(G4double sigma)
{
  if (sigma < 0.) {
    G4cerr << "BeamConfig: sigmaX must be >= 0. Keeping "
           << g_sigma_x / mm << " mm." << G4endl;
    return;
  }

  g_sigma_x = sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::SetSigmaY(G4double sigma)
{
  if (sigma < 0.) {
    G4cerr << "BeamConfig: sigmaY must be >= 0. Keeping "
           << g_sigma_y / mm << " mm." << G4endl;
    return;
  }

  g_sigma_y = sigma;
}

void BeamConfig::SetOffsetX(G4double offset) { g_offset_x = offset; }
void BeamConfig::SetOffsetY(G4double offset) { g_offset_y = offset; }
void BeamConfig::SetZ(G4double z) { g_z = z; }

void BeamConfig::SetProfileCommand(G4String profile) { SetProfile(profile); }
void BeamConfig::SetRadiusCommand(G4double radius) { SetRadius(radius); }
void BeamConfig::SetSigmaXCommand(G4double sigma) { SetSigmaX(sigma); }
void BeamConfig::SetSigmaYCommand(G4double sigma) { SetSigmaY(sigma); }
void BeamConfig::SetOffsetXCommand(G4double offset) { SetOffsetX(offset); }
void BeamConfig::SetOffsetYCommand(G4double offset) { SetOffsetY(offset); }
void BeamConfig::SetZCommand(G4double z) { SetZ(z); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::PrintConfigCommand()
{
  constexpr G4double fwhm_per_sigma = 2.3548200450309493;

  G4cout << "\n===== Beam runtime configuration =====" << G4endl
         << "  profile       = " << GetProfileName() << G4endl
         << "  radius        = " << GetRadius() / mm << " mm"
         << " (uniform-disk diameter = " << 2. * GetRadius() / mm << " mm)"
         << G4endl
         << "  sigmaX        = " << GetSigmaX() / mm << " mm"
         << " (FWHM = " << fwhm_per_sigma * GetSigmaX() / mm << " mm)"
         << G4endl
         << "  sigmaY        = " << GetSigmaY() / mm << " mm"
         << " (FWHM = " << fwhm_per_sigma * GetSigmaY() / mm << " mm)"
         << G4endl
         << "  offsetX       = " << GetOffsetX() / mm << " mm" << G4endl
         << "  offsetY       = " << GetOffsetY() / mm << " mm" << G4endl
         << "  source z      = " << GetZ() / mm << " mm" << G4endl
         << "  divergence    = 0 deg (beam direction remains +z)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void BeamConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(
      this,
      "/beam/",
      "Normal proton-beam profile controls");

  messenger->DeclareMethod(
      "profile",
      &BeamConfig::SetProfileCommand,
      "Transverse profile: point, uniformDisk or gaussian");

  messenger->DeclareMethodWithUnit(
      "radius",
      "mm",
      &BeamConfig::SetRadiusCommand,
      "Uniform-disk radius; radius = 0 is equivalent to a point beam");

  messenger->DeclareMethodWithUnit(
      "sigmaX",
      "mm",
      &BeamConfig::SetSigmaXCommand,
      "Gaussian standard deviation along global x");

  messenger->DeclareMethodWithUnit(
      "sigmaY",
      "mm",
      &BeamConfig::SetSigmaYCommand,
      "Gaussian standard deviation along global y");

  messenger->DeclareMethodWithUnit(
      "offsetX",
      "mm",
      &BeamConfig::SetOffsetXCommand,
      "Beam-centre offset along global x");

  messenger->DeclareMethodWithUnit(
      "offsetY",
      "mm",
      &BeamConfig::SetOffsetYCommand,
      "Beam-centre offset along global y");

  messenger->DeclareMethodWithUnit(
      "z",
      "mm",
      &BeamConfig::SetZCommand,
      "Primary proton source plane in global z");

  messenger->DeclareMethod(
      "printConfig",
      &BeamConfig::PrintConfigCommand,
      "Print the current beam-profile configuration");
}
