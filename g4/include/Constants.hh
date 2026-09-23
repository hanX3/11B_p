#ifndef Constants_H
#define Constants_H 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"

enum class H11BGammaResonance
{
  None = 0,
  Resonance165 = 165,
  Resonance675 = 675
};

enum class H11BGammaBranch
{
  None = 0,
  Gamma165ToGroundState = 1,
  Gamma165ToFirstExcitedState = 2,
  Gamma675ToGroundState = 3,
  Gamma675To4439State = 4,
  Gamma675To7654State = 5,
  Gamma675To12710State = 6,
  Gamma675To15110State = 7
};


enum class H11BGammaAngularMode
{
  None = 0,
  Isotropic = 1,
  FixedA1A2 = 2
};

// 675-keV alpha-decay final-state generator.
// LegacyLegendreA2A4 keeps the pre-existing sequential alpha+8Be(2+) model
// with an empirical A2/A4 secondary angular correlation.
// SymmetrizedCoherentL1L3 uses a Kuhlwein-2022-inspired symmetrized,
// coherent L=1/L=3 sequential 3-alpha amplitude for the 16.62-MeV 2- state.
enum class H11B675AlphaDecayModel
{
  LegacyLegendreA2A4 = 0,
  SymmetrizedCoherentL1L3 = 1
};


enum class DetectorGeometryMode
{
  CurrentChamberPortAligned = 0,
  IdealNoChamber = 1,
  ChamberWithExtraGammaWindows = 2
};

constexpr DetectorGeometryMode DetectorGeometryModeDefault = DetectorGeometryMode::CurrentChamberPortAligned;

// #define ReactionCout

// ROOT output is controlled at runtime through /output/ commands.
// Defaults are defined in OutputConfig.cc.

constexpr G4double WorldSizeX = 2. * m;
constexpr G4double WorldSizeY = 2. * m;
constexpr G4double WorldSizeZ = 2. * m;

// Chamber
constexpr G4double ChamberX = 257. * mm;
constexpr G4double ChamberY = 257. * mm;
constexpr G4double ChamberZ = 657. * mm;
constexpr G4double ChamberThickness = 10. * mm;

constexpr G4double FlangeXR = 104. * mm - ChamberThickness / 2.;
constexpr G4double FlangeYR = 77. * mm - ChamberThickness / 2.;
constexpr G4double FlangeZR = 52. * mm - ChamberThickness / 2.;

constexpr G4double FlangeXR2 = 126.5 * mm;
constexpr G4double FlangeYR2 = 101. * mm;
constexpr G4double FlangeZR2 = 60. * mm;

constexpr G4double FlangeXRH = 24.5 * mm;
constexpr G4double FlangeYRH = 24.5 * mm;
constexpr G4double FlangeZRH = 24.5 * mm;

// Target
constexpr G4double TargetR = 20. * mm;
constexpr G4double TargetThickness = 5. * mm;
constexpr G4double TargetZPos = 170. * mm;
const G4String TargetMaterial = "Natured_11B_low_density";

// TargetBacking
constexpr G4bool TargetBackingFlag = true;
constexpr G4double TargetBackingR = 20. * mm;
constexpr G4double TargetBackingThickness = 5. * mm;

// Si array
constexpr G4double SiAlShellThickness = 0.5 * mm;
constexpr G4double SiEnergyResolution = 0.02;
constexpr G4double SiEnergyThreshold = 20. * keV;

// LaBr3 array
constexpr G4double LaBr3AlShellThickness = 0.5 * mm;
constexpr G4double LaBr3EnergyResolution = 0.02;
constexpr G4double LaBr3EnergyThreshold = 20. * keV;

// HPGe array
constexpr G4double HPGeAlShellThickness = 0.5 * mm;
constexpr G4double HPGeEnergyResolution = 0.02;
constexpr G4double HPGeEnergyThreshold = 20. * keV;

// Beam
constexpr G4double BeamR = 5. * mm;
constexpr G4double BeamZ = 0. * mm;
constexpr G4double BeamEnergy = 165. * keV;

// Reaction
// 8Be(2+) line-shape parameters.
// Ex8Be2Plus is relative to the 8Be ground state and is added to the 8Be mass.
// Eaa8Be2Plus is relative to the 2-alpha threshold and is used in the
// energy-dependent width of the 8Be(2+) -> alpha + alpha decay.
constexpr G4double Ex8Be2Plus = 3.03 * MeV;
constexpr G4double Gamma8Be2Plus = 1.513 * MeV;
constexpr G4double Ex8BeGroundAbove2Alpha = 91.84 * keV;
constexpr G4double Eaa8Be2Plus = Ex8Be2Plus + Ex8BeGroundAbove2Alpha;
constexpr G4int L8Be2PlusAlphaAlpha = 2;

// Backward-compatible names used by older code sections.
constexpr G4double Ex8Be = Ex8Be2Plus;
constexpr G4double Ex8BeGamma = Gamma8Be2Plus;

// p + 11B Breit-Wigner parameters from TUNL Table 12.20.
// Energies used by the cross-section code are in the center-of-mass system.
// The TUNL table lists proton beam energies Ep in the lab system; therefore
// Ecm = 11/12 * Ep is used for a proton incident on a stationary 11B target.
// Center-only scan against the Tentori first peak gives Er = 147.95 keV
// in the center-of-mass system, aligning the generated BW maximum with the
// evaluated total-cross-section peak near Ep(lab) = 161.55 keV.
constexpr G4double H11B165ResonanceEnergy = 147.95 * keV; // Ex ~= 16.106 MeV
constexpr G4double H11B165TotalWidth = 5.3 * keV;
constexpr G4double H11B165ProtonWidth = 0.0215 * keV;
constexpr G4double H11B165Alpha0Width = 0.26 * keV;
constexpr G4double H11B165Alpha1Width = 5.0 * keV;
constexpr G4double H11B165Gamma0Width = 0.66 * eV;
constexpr G4double H11B165Gamma1Width = 18.0 * eV;
constexpr G4double H11B165GammaWidth = H11B165Gamma0Width + H11B165Gamma1Width;
constexpr G4double H11B165SpinStatFactor = 5.0 / 8.0;

// The 16.576-MeV 2- resonance cannot decay through the alpha + 8Be(g.s.)
// alpha0 channel by parity conservation, so the 675-keV alpha0 width is fixed
// to zero in this model.
constexpr G4double H11B675Alpha0Width = 0.0 * keV;
constexpr G4double H11B675Alpha1Width = 150.0 * keV;

// The entrance-channel proton width is always calculated using the full
// Coulomb penetrability ratio P_l(E) / P_l(E_r).

// Entrance-channel orbital angular momentum for p + 11B.
constexpr G4int H11B165EntranceOrbitalL = 1; // likely p-wave

// Use full Coulomb penetrability for 8Be(2+) -> alpha + alpha line shape.
constexpr G4bool H11BUseFullCoulombPenetrabilityFor8Be2Plus = true;

// Use full Coulomb penetrability for the first breakup
//   12C* -> alpha + 8Be(2+)
// when sampling the 8Be(2+) excitation energy in the alpha1 branch.
// This is separate from the 8Be(2+) -> alpha + alpha penetrability above.
constexpr G4bool H11BUseFullCoulombPenetrabilityForAlpha8BeFirstBreakup = true;

// Primary-alpha angular distribution for the 165-keV resonance first breakup
//   12C*(16.11, 2+) -> alpha + 8Be.
// Becker 1987 fits the low-energy alpha0/alpha1 primary angular distributions
// with
//   W(theta) = 1 + a1*P1(cos(theta)) + a2*P2(cos(theta)).
// The code convention is theta relative to the incident proton beam axis in
// the 12C center-of-mass frame. If Becker coefficients digitized from the
// figures are defined relative to a reversed 11B beam axis, odd coefficients
// must be sign-flipped before use: a1_code=-a1_becker, a2_code=a2_becker.
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// IMPORTANT: the default remains isotropic until Becker Fig. 10/11 coefficients
// are digitized and validated. When enabled from macro, the same A1/A2
// coefficients are used for the sampled 165-keV primary alpha direction.
constexpr G4bool H11BDefaultEnable165PrimaryAngularDistribution = false;
constexpr G4double H11BDefault165PrimaryAngularA1 = 0.0;
constexpr G4double H11BDefault165PrimaryAngularA2 = 0.0;

// Exit-channel angular-correlation defaults.
// These affect only the 8Be(2+) sequential branch, i.e. the alpha1 channel.
// Treado 1972 gives the 8Be(2+) breakup correlation as
//   W(chi) = 1 + A2*P2(cos chi) + A4*P4(cos chi)
// when the small fitted symmetry-axis shift is ignored.
// The isotropic modes remain available in macros for diagnostics.
constexpr G4bool H11BDefaultEnable165Alpha1SecondaryAngularCorrelation = true;
constexpr G4double H11BDefault165Alpha1SecondaryA2 = -0.489;
constexpr G4double H11BDefault165Alpha1SecondaryA4 = 0.647;
constexpr G4bool H11BDefaultEnable675PrimaryAngularDistribution = false;
constexpr G4double H11BDefault675PrimaryAngularA1 = 0.0;
constexpr G4double H11BDefault675PrimaryAngularA2 = 0.0;
constexpr G4bool H11BDefaultEnable675Alpha1SecondaryAngularCorrelation = true;
constexpr G4double H11BDefault675Alpha1SecondaryA2 = -1.064;
constexpr G4double H11BDefault675Alpha1SecondaryA4 = 0.180;

// 165-keV resonance: 12C*(16.11, 2+) -> alpha + 8Be(2+).
// The alpha1 first-breakup penetrability is treated as d-wave dominated: L = 2.
constexpr G4int H11B165Alpha1ExitOrbitalL = 2;

// 675-keV resonance: 12C*(16.57/16.62, 2-) -> alpha + 8Be(2+).
// The primary-alpha direction can optionally be sampled with the same
// A1/A2 Legendre form used for the 165-keV primary alpha.  The internal
// 8Be(2+) -> alpha + alpha correlation can optionally use the Treado-1972
// A2/A4 coefficients.  The L=1/L=3 constants below are retained only for
// first-breakup penetrability weighting and future model tests.
constexpr G4int H11B675Alpha1ExitOrbitalL1 = 1;
constexpr G4int H11B675Alpha1ExitOrbitalL3 = 3;
constexpr G4double H11B675ExitL1Fraction = 0.76;
constexpr G4double H11B675ExitL13Phase = 0.67 * 2.0 * 3.14159265358979323846;
constexpr G4bool H11B675UseCoherentL13Interference = false;

// Strict 675-keV alpha-decay defaults.  The model follows the 2022 PLB
// exclusive-decay analysis at the event-generator level: alpha+8Be(2+)
// sequential decay, coherent L=1/L=3 primary-emission amplitudes, and
// symmetrization over the three identical alpha particles.
constexpr H11B675AlphaDecayModel H11BDefault675AlphaDecayModel = H11B675AlphaDecayModel::SymmetrizedCoherentL1L3;
constexpr G4double H11B675StrictDefaultL1Fraction = 0.76;
constexpr G4double H11B675StrictDefaultL13Phase = 0.67 * 2.0 * 3.14159265358979323846;
constexpr G4bool H11B675StrictDefaultCoherentL13 = true;
constexpr G4bool H11B675StrictDefaultPermutationSymmetrized = true;
constexpr G4double H11B675StrictDefault8BeLambdaEnergy = 3037.0 * keV;
constexpr G4double H11B675StrictDefault8BeReducedWidthSquared = 1075.0 * keV;
constexpr G4double H11B675StrictDefaultWeightMaxSafetyFactor = 3.0;
constexpr G4int H11B675StrictDefaultWeightMaxScanCandidates = 2000;
constexpr G4int H11B675StrictDefaultMaxSamplingAttempts = 10000;

// 12C gamma capture channels following 11B(p,gamma)12C.
// Gamma capture is an independent exit channel competing with 3alpha; it is
// not included in the phenomenological directdecay yield.
constexpr G4bool H11BDefault165GammaCaptureEnabled = true;
constexpr G4double H11BDefaultGammaBiasFactor = 1.0;

// Optional fallback scales relative to the 165 alpha model cross section.
// The production cross section uses the literature-based gamma partial widths
// above; these constants are retained only for sensitivity studies.
constexpr G4double H11B165Gamma0ToAlphaScale = 1.2e-4;
constexpr G4double H11B165Gamma1ToAlphaScale = 3.4e-3;

// 165 gamma0 angular distribution in fixed Legendre A1/A2 form:
//   W(theta) = 1 + a1*P1(cos(theta)) + a2*P2(cos(theta)).
// The default values below are equivalent, up to normalization, to the
// Craig-1956 form W(theta)=1-0.19*cos(theta)+0.21*cos(theta)^2.
constexpr G4bool H11BDefaultEnable165Gamma0AngularDistribution = true;
constexpr G4double H11B165Gamma0AngularA1Default = -0.17757009345794392;
constexpr G4double H11B165Gamma0AngularA2Default = 0.13084112149532712;

// 675-keV primary gamma angular distribution.  TUNL/Ajzenberg-Selove
// describes the 675-keV resonance gamma1 exit channel as near isotropic;
// therefore the default remains isotropic unless explicit A1/A2 coefficients
// are supplied for a sensitivity test.
constexpr G4bool H11BDefaultEnable675GammaAngularDistribution = false;
constexpr G4double H11B675GammaAngularA1Default = 0.0;
constexpr G4double H11B675GammaAngularA2Default = 0.0;

// 675-keV gamma is a phenomenological sensitivity model.
// 675-keV gamma branches are sampled internally from the relative line table.
constexpr G4bool H11BDefault675GammaCaptureEnabled = true;
constexpr G4double H11B675GammaToAlphaScale = 1.0e-5;
constexpr G4bool H11B675IncludeUpperLimitGammaLines = false;

// 12C level energies used by gamma cascade kinematics.
constexpr G4double C12LevelGround = 0.0 * MeV;
constexpr G4double C12Level4439 = 4.439 * MeV;
constexpr G4double C12Level7654 = 7.654 * MeV;
constexpr G4double C12Level9641 = 9.641 * MeV;
constexpr G4double C12Level12710 = 12.710 * MeV;
constexpr G4double C12Level15110 = 15.110 * MeV;
constexpr G4double C12Resonance165Ex = 16.106 * MeV;
constexpr G4double C12Resonance675Ex = 16.57 * MeV;

// Backward-compatible flag name.  Gamma final states are generated explicitly;
// alpha-model Breit-Wigner exit widths must not include gamma widths.
constexpr G4bool H11BIncludeGammaChannelInCrossSection = false;

constexpr G4double StepMax4Proton = 2 * nm;

#endif
