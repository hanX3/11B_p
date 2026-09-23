#ifndef SiArrayConfig_H
#define SiArrayConfig_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <memory>

class G4GenericMessenger;

// Runtime configuration for the silicon array.
//
// Default geometry in this version:
//   * 12 W1 DSSDs form the equatorial drum.
//   * Twelve wedge DSSSDs form a forward CAKE-like lampshade.
//   * Twelve mirrored wedge DSSSDs form a backward lampshade.
//   * One S3 is placed behind the central aperture of each lampshade.
//
// The lampshade plane is defined directly by its centre distance and centre
// polar angle.  Its long local axis is tangent to a sphere centred on the
// target, producing a visibly inclined lampshade rather than a nearly flat
// annulus.
class SiArrayConfig
{
public:
  SiArrayConfig();
  ~SiArrayConfig();

  // ---- W1 drum readout ----
  static G4int GetBarrelStripsZ();
  static G4int GetBarrelStripsPhi();

  // ---- S3 readout ----
  static G4int GetAnnularRings();
  static G4int GetAnnularSectors();

  // ---- Lampshade wedge-DSSD readout ----
  static G4int GetLampshadeRings();
  static G4int GetLampshadeSectors();

  // ---- S3 geometry ----
  static G4bool GetEnableForwardAnnular();
  static G4bool GetEnableBackwardAnnular();
  static G4double GetForwardDistance();
  static G4double GetBackwardDistance();

  // ---- W1 drum geometry ----
  static G4int GetDrumModules();
  static G4double GetDrumInscribedRadius();
  static G4bool GetEnableDrum();

  // ---- CAKE-like lampshade geometry ----
  static G4bool GetEnableForwardCap();
  static G4bool GetEnableBackwardCap();
  static G4int GetLampshadeModules();
  static G4double GetLampshadeCenterDistance();
  static G4double GetLampshadeCenterAngle();
  static G4double GetLampshadeLength();
  static G4double GetLampshadeModuleSpan();
  static G4double GetLampshadeThickness();

  // Derived narrow-edge and wide-edge geometry. Distances are positive axial
  // distances from the target centre for one end; the backward end is mirrored.
  static G4double GetLampshadeInnerRadius();
  static G4double GetLampshadeInnerDistance();
  static G4double GetLampshadeOuterRadius();
  static G4double GetLampshadeOuterDistance();
  static G4double GetLampshadeInnerAngle();
  static G4double GetLampshadeOuterAngle();
  static G4double GetLampshadeInnerHalfWidth();
  static G4double GetLampshadeOuterHalfWidth();

  // Historical aliases retained to avoid breaking older source code/macros.
  static G4double GetLampshadeHoleApothem();
  static G4double GetCapFoldAngle();
  static G4double GetCapClearance();
  static G4double GetCapStagger();

  static void SetBarrelStripsZ(G4int n);
  static void SetBarrelStripsPhi(G4int n);
  static void SetAnnularRings(G4int n);
  static void SetAnnularSectors(G4int n);
  static void SetLampshadeRings(G4int n);
  static void SetLampshadeSectors(G4int n);

  static void SetEnableForwardAnnular(G4bool enabled);
  static void SetEnableBackwardAnnular(G4bool enabled);
  static void SetForwardDistance(G4double d);
  static void SetBackwardDistance(G4double d);

  static void SetDrumModules(G4int n);
  static void SetEnableDrum(G4bool enabled);

  static void SetEnableForwardCap(G4bool enabled);
  static void SetEnableBackwardCap(G4bool enabled);
  static void SetLampshadeModules(G4int n);
  static void SetLampshadeCenterDistance(G4double d);
  static void SetLampshadeCenterAngle(G4double angle);
  static void SetLampshadeLength(G4double d);
  static void SetLampshadeModuleSpan(G4double angle);
  static void SetLampshadeThickness(G4double d);

  // Obsolete geometry setters retained for macro compatibility.
  static void SetLampshadeHoleApothem(G4double d);
  static void SetLampshadeInnerDistance(G4double d);
  static void SetLampshadeOuterAngle(G4double angle);
  static void SetCapFoldAngle(G4double angle);
  static void SetCapClearance(G4double d);
  static void SetCapStagger(G4double d);

  // ---- Detector-test primary modes (default OFF) ----
  static G4bool GetForwardAlphaTestBeam();
  static G4double GetForwardAlphaTestBeamEnergy();
  static void SetForwardAlphaTestBeam(G4bool enabled);
  static void SetForwardAlphaTestBeamEnergy(G4double energy);

  static G4bool GetFixedAlphaTestBeam();
  static G4double GetFixedAlphaTestBeamEnergy();
  static G4double GetFixedAlphaTestBeamTheta();
  static G4double GetFixedAlphaTestBeamPhi();
  static void SetFixedAlphaTestBeam(G4bool enabled);
  static void SetFixedAlphaTestBeamEnergy(G4double energy);
  static void SetFixedAlphaTestBeamTheta(G4double theta);
  static void SetFixedAlphaTestBeamPhi(G4double phi);

  // ---- Strip IDs from local hit coordinates ----
  static G4int BarrelPhiStripId(G4double local_x, G4double half_x);
  static G4int BarrelZStripId(G4double local_y, G4double half_y);
  static G4int AnnularSectorStripId(G4double local_x, G4double local_y);
  static G4int AnnularRingStripId(G4double local_x,
                                 G4double local_y,
                                 G4double r_inner,
                                 G4double r_outer);
  static G4double LampshadeHalfWidthAtY(G4double local_y,
                                       G4double half_length,
                                       G4double inner_half_width,
                                       G4double outer_half_width);
  static G4int LampshadeSectorStripId(G4double local_x,
                                     G4double local_y,
                                     G4double half_length,
                                     G4double inner_half_width,
                                     G4double outer_half_width);
  static G4int LampshadeRingStripId(G4double local_y, G4double half_length);

  static G4int BarrelReadoutChannelCount();
  static G4int AnnularReadoutChannelCount();
  static G4int LampshadeReadoutChannelCount();

  static void SelfTestCopyNoEncoding();

  // Messenger wrappers.
  void SetBarrelStripsZCmd(G4int n);
  void SetBarrelStripsPhiCmd(G4int n);
  void SetAnnularRingsCmd(G4int n);
  void SetAnnularSectorsCmd(G4int n);
  void SetLampshadeRingsCmd(G4int n);
  void SetLampshadeSectorsCmd(G4int n);
  void SetEnableForwardAnnularCmd(G4bool enabled);
  void SetEnableBackwardAnnularCmd(G4bool enabled);
  void SetForwardDistanceCmd(G4double d);
  void SetBackwardDistanceCmd(G4double d);
  void SetDrumModulesCmd(G4int n);
  void SetEnableDrumCmd(G4bool enabled);
  void SetEnableForwardCapCmd(G4bool enabled);
  void SetEnableBackwardCapCmd(G4bool enabled);
  void SetLampshadeModulesCmd(G4int n);
  void SetLampshadeCenterDistanceCmd(G4double d);
  void SetLampshadeCenterAngleCmd(G4double angle);
  void SetLampshadeLengthCmd(G4double d);
  void SetLampshadeModuleSpanCmd(G4double angle);
  void SetLampshadeThicknessCmd(G4double d);
  void SetLampshadeHoleApothemCmd(G4double d);
  void SetLampshadeInnerDistanceCmd(G4double d);
  void SetLampshadeOuterAngleCmd(G4double angle);
  void SetCapFoldAngleCmd(G4double angle);
  void SetCapClearanceCmd(G4double d);
  void SetCapStaggerCmd(G4double d);
  void SetForwardAlphaTestBeamCmd(G4bool enabled);
  void SetForwardAlphaTestBeamEnergyCmd(G4double energy);
  void SetFixedAlphaTestBeamCmd(G4bool enabled);
  void SetFixedAlphaTestBeamEnergyCmd(G4double energy);
  void SetFixedAlphaTestBeamThetaCmd(G4double theta);
  void SetFixedAlphaTestBeamPhiCmd(G4double phi);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
