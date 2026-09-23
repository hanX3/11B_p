#include "VirtualSphereConfig.hh"

#include "Constants.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

namespace {
G4bool g_enabled = true;
G4double g_radius = 25. * mm;
G4double g_thickness = 1. * um;
G4bool g_save_electrons = false;
G4bool g_save_optical_photons = false;
G4double g_min_kinetic_energy = 0.;
} // namespace

VirtualSphereConfig::VirtualSphereConfig()
{
  DefineCommands();
}

VirtualSphereConfig::~VirtualSphereConfig() = default;

G4bool VirtualSphereConfig::GetEnabled() { return g_enabled; }
G4double VirtualSphereConfig::GetRadius() { return g_radius; }
G4double VirtualSphereConfig::GetThickness() { return g_thickness; }
G4bool VirtualSphereConfig::GetSaveElectrons() { return g_save_electrons; }
G4bool VirtualSphereConfig::GetSaveOpticalPhotons() { return g_save_optical_photons; }
G4double VirtualSphereConfig::GetMinKineticEnergy() { return g_min_kinetic_energy; }

void VirtualSphereConfig::SetEnabled(G4bool enabled) { g_enabled = enabled; }

void VirtualSphereConfig::SetRadius(G4double radius)
{
  if (radius > TargetR) g_radius = radius;
}

void VirtualSphereConfig::SetThickness(G4double thickness)
{
  if (thickness > 0.) g_thickness = thickness;
}

void VirtualSphereConfig::SetSaveElectrons(G4bool enabled) { g_save_electrons = enabled; }
void VirtualSphereConfig::SetSaveOpticalPhotons(G4bool enabled) { g_save_optical_photons = enabled; }

void VirtualSphereConfig::SetMinKineticEnergy(G4double energy)
{
  g_min_kinetic_energy = energy >= 0. ? energy : 0.;
}

void VirtualSphereConfig::SetEnabledCommand(G4bool enabled) { SetEnabled(enabled); }
void VirtualSphereConfig::SetRadiusCommand(G4double radius) { SetRadius(radius); }
void VirtualSphereConfig::SetThicknessCommand(G4double thickness) { SetThickness(thickness); }
void VirtualSphereConfig::SetSaveElectronsCommand(G4bool enabled) { SetSaveElectrons(enabled); }
void VirtualSphereConfig::SetSaveOpticalPhotonsCommand(G4bool enabled) { SetSaveOpticalPhotons(enabled); }
void VirtualSphereConfig::SetMinKineticEnergyCommand(G4double energy) { SetMinKineticEnergy(energy); }

void VirtualSphereConfig::PrintConfigCommand()
{
  G4cout << "\n===== Virtual sphere configuration =====" << G4endl
         << "  enabled             = " << (GetEnabled() ? "true" : "false") << G4endl
         << "  radius              = " << GetRadius() / mm << " mm" << G4endl
         << "  thickness           = " << GetThickness() / um << " um" << G4endl
         << "  saveElectrons       = " << (GetSaveElectrons() ? "true" : "false") << G4endl
         << "  saveOpticalPhotons  = " << (GetSaveOpticalPhotons() ? "true" : "false") << G4endl
         << "  minKineticEnergy    = " << GetMinKineticEnergy() / keV << " keV" << G4endl;
}

void VirtualSphereConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/virtualSphere/",
                                                    "Virtual target-surrounding scoring shell");

  messenger->DeclareMethod("enabled", &VirtualSphereConfig::SetEnabledCommand,
                           "Enable/disable the virtual scoring sphere (set before /run/initialize)");
  messenger->DeclareMethodWithUnit("radius", "mm", &VirtualSphereConfig::SetRadiusCommand,
                                   "Inner radius; must exceed the target radius");
  messenger->DeclareMethodWithUnit("thickness", "um", &VirtualSphereConfig::SetThicknessCommand,
                                   "Numerical shell thickness; material is chamber vacuum");
  messenger->DeclareMethod("saveElectrons", &VirtualSphereConfig::SetSaveElectronsCommand,
                           "Record e-/e+ crossings");
  messenger->DeclareMethod("saveOpticalPhotons", &VirtualSphereConfig::SetSaveOpticalPhotonsCommand,
                           "Record optical-photon crossings");
  messenger->DeclareMethodWithUnit("minKineticEnergy", "keV",
                                   &VirtualSphereConfig::SetMinKineticEnergyCommand,
                                   "Minimum kinetic energy for a saved crossing");
  messenger->DeclareMethod("printConfig", &VirtualSphereConfig::PrintConfigCommand,
                           "Print current virtual-sphere configuration");
}
