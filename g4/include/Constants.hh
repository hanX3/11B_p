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
constexpr G4int MASK = 0b011;

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
// 8Be(2+) line-shape parameters.
// Ex8Be2Plus is relative to the 8Be ground state and is added to the 8Be mass.
// Eaa8Be2Plus is relative to the 2-alpha threshold and is used in the
// energy-dependent width of the 8Be(2+) -> alpha + alpha decay.
constexpr G4double Ex8Be2Plus = 3.03 *MeV;
constexpr G4double Gamma8Be2Plus = 1.513 *MeV;
constexpr G4double Ex8BeGroundAbove2Alpha = 91.84 *keV;
constexpr G4double Eaa8Be2Plus = Ex8Be2Plus + Ex8BeGroundAbove2Alpha;
constexpr G4int L8Be2PlusAlphaAlpha = 2;

// Backward-compatible names used by older code sections.
constexpr G4double Ex8Be = Ex8Be2Plus;
constexpr G4double Ex8BeGamma = Gamma8Be2Plus;

// p + 11B Breit-Wigner parameters from TUNL Table 12.20.
// Energies used by the cross-section code are in the center-of-mass system.
// The TUNL table lists proton beam energies Ep in the lab system; therefore
// Ecm = 11/12 * Ep is used for a proton incident on a stationary 11B target.
constexpr G4double H11B165ResonanceEnergy = 148.3 *keV;  // TUNL: Ep = 0.162 MeV, Ex = 16.106 MeV
constexpr G4double H11B165TotalWidth = 5.3 *keV;
constexpr G4double H11B165ProtonWidth = 0.0215 *keV;
constexpr G4double H11B165Alpha0Width = 0.26 *keV;
constexpr G4double H11B165Alpha1Width = 5.0 *keV;
constexpr G4double H11B165Gamma0Width = 0.59 *eV;
constexpr G4double H11B165Gamma1Width = 12.8 *eV;
constexpr G4double H11B165GammaWidth = H11B165Gamma0Width + H11B165Gamma1Width;
constexpr G4double H11B165SpinStatFactor = 5.0/8.0;

constexpr G4double H11B675ResonanceEnergy = 618.75 *keV;  // TUNL: Ep = 0.675 MeV, Ex = 16.576 MeV
constexpr G4double H11B675TotalWidth = 300.0 *keV;
constexpr G4double H11B675ProtonWidth = 150.0 *keV;
// TUNL gives Gamma_alpha0 < 0.27 keV. For a pure 16.576-MeV 2- resonance,
// alpha0 is parity-forbidden in the alpha + 8Be(g.s.) channel, so the default
// generated resonant alpha0 branch is set to zero. Change the flag below to
// true for an upper-limit sensitivity test.
constexpr G4bool H11B675UseAlpha0UpperLimit = false;
constexpr G4double H11B675Alpha0UpperLimit = 0.27 *keV;
constexpr G4double H11B675Alpha0Width = H11B675UseAlpha0UpperLimit ? H11B675Alpha0UpperLimit : 0.0 *keV;
constexpr G4double H11B675Alpha1Width = 150.0 *keV;
constexpr G4double H11B675Gamma0UpperLimit = 0.4 *eV;
constexpr G4double H11B675Gamma0Width = H11B675Gamma0UpperLimit;
constexpr G4double H11B675Gamma1Width = 8.0 *eV;
constexpr G4double H11B675GammaWidth = H11B675Gamma0Width + H11B675Gamma1Width;
constexpr G4double H11B675SpinStatFactor = 5.0/8.0;

// Use full Coulomb penetrability instead of simple power-law approximation.
constexpr G4bool H11BUseFullCoulombPenetrabilityInCrossSection = true;

// Entrance-channel orbital angular momentum for p + 11B.
constexpr G4int H11B165EntranceOrbitalL = 1;  // likely p-wave
constexpr G4int H11B675EntranceOrbitalL = 0;  // likely s-wave, check with adopted resonance assignment

// Use full Coulomb penetrability for 8Be(2+) -> alpha + alpha line shape.
constexpr G4bool H11BUseFullCoulombPenetrabilityFor8Be2Plus = true;

// Primary-alpha angular distribution for the first breakup
//   12C*(16.11, 2+) -> alpha + 8Be.
// The implemented form is W(theta) = 1 + a2*P2(cos(theta)) + a4*P4(cos(theta)),
// where theta is measured relative to the incident-proton/beam axis in the
// 12C center-of-mass frame.  The coefficients are branch dependent because
// alpha0 and alpha1 can have different experimental angular distributions.
//
// IMPORTANT: these coefficients are not fixed by L_out=2 alone; they depend on
// entrance-channel spin amplitudes / alignment.  Set them from data or a more
// complete R-matrix/alignment calculation.  The default zeros preserve the old
// isotropic primary-alpha distribution until physical coefficients are supplied.
constexpr G4bool H11B165UsePrimaryAngularDistribution = false;
constexpr G4double H11B165Alpha0PrimaryA2 = 0.0;
constexpr G4double H11B165Alpha0PrimaryA4 = 0.0;
constexpr G4double H11B165Alpha1PrimaryA2 = 0.0;
constexpr G4double H11B165Alpha1PrimaryA4 = 0.0;

// Exit-channel angular-correlation switches.
// These affect only the 8Be(2+) sequential branch, i.e. alpha1 channel.
constexpr G4bool H11B165UseExitAngularCorrelation = false;
constexpr G4bool H11B675UseExitAngularCorrelation = false;

// 165-keV resonance: 12C*(16.11, 2+) -> alpha + 8Be(2+).
// The alpha1 branch is treated as d-wave dominated: L = 2.
constexpr G4int H11B165Alpha1ExitOrbitalL = 2;

// 675-keV resonance: 12C*(16.57/16.62, 2-) -> alpha + 8Be(2+).
// Coherent L=1/L=3 mixture.  Set L1Fraction=1 for pure L=1 and 0 for pure L=3.
constexpr G4int H11B675Alpha1ExitOrbitalL1 = 1;
constexpr G4int H11B675Alpha1ExitOrbitalL3 = 3;
constexpr G4double H11B675ExitL1Fraction = 0.76;
constexpr G4double H11B675ExitL13Phase = 0.67 * 2.0 * 3.14159265358979323846;
constexpr G4bool H11B675UseCoherentL13Interference = false;


// The current H11BReaction final state still generates 3 alpha particles.
// Keep this flag false unless a gamma final-state generator is also added.
constexpr G4bool H11BIncludeGammaChannelInCrossSection = false;

constexpr G4double StepMax4Proton = 40 *nm;


#endif
