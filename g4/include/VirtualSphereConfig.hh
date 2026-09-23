#ifndef VirtualSphereConfig_H
#define VirtualSphereConfig_H 1

#include "globals.hh"

#include <memory>

class G4GenericMessenger;

// Runtime controls for the non-physical vacuum scoring shell surrounding the
// target. Geometry-affecting commands must be issued before /run/initialize.
class VirtualSphereConfig
{
public:
  VirtualSphereConfig();
  ~VirtualSphereConfig();

  static G4bool GetEnabled();
  static G4double GetRadius();
  static G4double GetThickness();
  static G4bool GetSaveElectrons();
  static G4bool GetSaveOpticalPhotons();
  static G4double GetMinKineticEnergy();

  static void SetEnabled(G4bool enabled);
  static void SetRadius(G4double radius);
  static void SetThickness(G4double thickness);
  static void SetSaveElectrons(G4bool enabled);
  static void SetSaveOpticalPhotons(G4bool enabled);
  static void SetMinKineticEnergy(G4double energy);

  void SetEnabledCommand(G4bool enabled);
  void SetRadiusCommand(G4double radius);
  void SetThicknessCommand(G4double thickness);
  void SetSaveElectronsCommand(G4bool enabled);
  void SetSaveOpticalPhotonsCommand(G4bool enabled);
  void SetMinKineticEnergyCommand(G4double energy);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
