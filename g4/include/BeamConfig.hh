#ifndef BeamConfig_H
#define BeamConfig_H 1

#include "globals.hh"

#include <memory>

class G4GenericMessenger;

enum class BeamProfile
{
  Point,
  UniformDisk,
  Gaussian
};

// Runtime controls for the normal proton-beam transverse profile.
//
// These settings do not affect the fixed-alpha or forward-alpha detector-test
// modes controlled under /si/.  The defaults reproduce the previous behavior:
// a uniform disk of radius BeamR centred at (0, 0, BeamZ).
class BeamConfig
{
public:
  BeamConfig();
  ~BeamConfig();

  static BeamProfile GetProfile();
  static G4String GetProfileName();

  static G4double GetRadius();
  static G4double GetSigmaX();
  static G4double GetSigmaY();
  static G4double GetOffsetX();
  static G4double GetOffsetY();
  static G4double GetZ();

  static void SetProfile(const G4String& profile);
  static void SetRadius(G4double radius);
  static void SetSigmaX(G4double sigma);
  static void SetSigmaY(G4double sigma);
  static void SetOffsetX(G4double offset);
  static void SetOffsetY(G4double offset);
  static void SetZ(G4double z);

  void SetProfileCommand(G4String profile);
  void SetRadiusCommand(G4double radius);
  void SetSigmaXCommand(G4double sigma);
  void SetSigmaYCommand(G4double sigma);
  void SetOffsetXCommand(G4double offset);
  void SetOffsetYCommand(G4double offset);
  void SetZCommand(G4double z);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
