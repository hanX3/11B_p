#ifndef Constants_h
#define Constants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

#define DATAPATH "../data"  // no "/" in th end

// #define ReactionCout

// 0b0001 event data
// 0b0010 track data
// 0b0100 step data
constexpr G4int MASK = 0b111;

constexpr G4double WorldSizeX = 2. *m;
constexpr G4double WorldSizeY = 2. *m;
constexpr G4double WorldSizeZ = 2. *m;

// Chamber
constexpr G4double ChamberX = 257. *mm;
constexpr G4double ChamberY = 257. *mm;
constexpr G4double ChamberZ = 657. *mm;
constexpr G4double ChamberThickness = 10. *mm;

constexpr G4double FlangeXR = 104. *mm - ChamberThickness/2.;
constexpr G4double FlangeYR = 77. *mm - ChamberThickness/2.;
constexpr G4double FlangeZR = 52. *mm - ChamberThickness/2.;

constexpr G4double FlangeXR2 = 126.5 *mm;
constexpr G4double FlangeYR2 = 101. *mm;
constexpr G4double FlangeZR2 = 60. *mm;

constexpr G4double FlangeXRH = 24.5 *mm;
constexpr G4double FlangeYRH = 24.5 *mm;
constexpr G4double FlangeZRH = 24.5 *mm;

// Target
constexpr G4double TargetR = 20. *mm;
constexpr G4double TargetThickness = 5. *mm;
constexpr G4double TargetZPos = 170. *mm;
const G4String TargetMaterial = "Natured_11B_low_density";

// TargetBacking
constexpr G4bool TargetBackingFlag = true;
constexpr G4double TargetBackingR = 20. *mm;
constexpr G4double TargetBackingThickness = 5. *mm;

// Si array
constexpr G4double SiAlShellThickness = 0.5 *mm;
constexpr G4double SiEnergyResolution = 0.02;
constexpr G4double SiEnergyThreshold = 20. *keV;

// LaBr3 array
constexpr G4double LaBr3AlShellThickness = 0.5 *mm;
constexpr G4double LaBr3EnergyResolution = 0.02;
constexpr G4double LaBr3EnergyThreshold = 20. *keV;

// HPGe array
constexpr G4double HPGeAlShellThickness = 0.5 *mm;
constexpr G4double HPGeEnergyResolution = 0.02;
constexpr G4double HPGeEnergyThreshold = 20. *keV;

// Beam
constexpr G4double BeamR = 5. *mm;
constexpr G4double BeamZ = 0. *mm;
constexpr G4double BeamEnergy = 165. *keV;

// Reaction
constexpr G4double Ex8Be = 3.03 *MeV;
constexpr G4double Ex8BeGamma = 1.513 *MeV;

constexpr G4double StepMax4Proton = 40 *nm;


#endif
