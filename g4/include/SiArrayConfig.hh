#ifndef SiArrayConfig_H
#define SiArrayConfig_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <memory>

class G4GenericMessenger;

// Centralized, runtime-configurable geometry / strip-segmentation parameters
// for the Si (DSSD-style) array.  Exposed under /si/ so that sensitivity scans
// can be run without recompiling.
//
// ---- "Connected drum" baseline (all-commercial W1/S3 parts) ----
//   * Si_Drum        (ring_id 1): 12 W1-type DSSDs (50x50 mm, 16x16 strips)
//     forming a regular dodecagonal prism around the target.  The inscribed
//     radius is DERIVED from the module count so the faces always close:
//         r_drum = 25 mm / tan(pi / n_modules)   (n=12 -> 93.30 mm)
//     Polar coverage theta ~ 75-105 deg.
//   * Si_ForwardCap  (ring_id 4) / Si_BackwardCap (ring_id 5): 12 W1 each,
//     hinged on the drum end edges and folded inward by capFoldAngle.
//     Square plates folded inward would physically intersect, so each cap is
//     built as two fish-scale sub-rings of 6: even module ids hinge on the
//     drum edge, odd module ids hinge capStagger further out in radius and
//     are rotated half a face period, covering the V-gaps of the inner
//     sub-ring.  Fold 45 deg -> theta ~ 44-75 (forward) / 105-136 (backward).
//   * Si_Forward/BackwardAnnular (ring_id 3/2): S3-type annular DSSDs
//     (active r = 11-35 mm, 24 rings x 32 sectors) at +/-36 mm,
//     theta ~ 17-44 / 136-163, sealing the ends.
//
// IMPORTANT timing note: strip counts are consumed inside SiSD::ProcessHits.
// A single sensitive silicon solid produces two independent strip-hit
// collections, one for each DSSD face, so strip counts take effect immediately
// and do NOT require a geometry rebuild.  All
// module counts, fold angle, stagger, annular distances and enable flags DO
// affect the built geometry; set them before /run/initialize, or set them and
// issue /run/reinitializeGeometry.
class SiArrayConfig
{
public:
  SiArrayConfig();
  ~SiArrayConfig();

  // ---- Strip counts (consumed at hit time, no rebuild needed) ----
  // W1 modules (drum + caps): "z" = polar-direction strips (local y),
  // "phi" = azimuthal-direction strips (local x).  16x16 for a real W1.
  static G4int GetBarrelStripsZ();
  static G4int GetBarrelStripsPhi();
  static G4int GetAnnularRings();     // radial rings on the S3 annulars (24)
  static G4int GetAnnularSectors();   // angular sectors on the S3 annulars (32)

  // ---- Geometry-affecting parameters (consumed at build time) ----
  static G4bool GetEnableForwardAnnular();
  static G4double GetForwardDistance();   // |z - target| of forward annular
  static G4double GetBackwardDistance();  // |z - target| of backward annular

  static G4int GetDrumModules();          // W1 faces of the dodecagonal drum
  static G4double GetDrumInscribedRadius(); // DERIVED: 25mm / tan(pi/n)
  static G4double GetCapFoldAngle();      // inward fold of the cap plates
  static G4double GetCapStagger();        // radial offset of the outer fish-scale sub-ring
  static G4bool GetEnableDrum();
  static G4bool GetEnableForwardCap();
  static G4bool GetEnableBackwardCap();

  static void SetBarrelStripsZ(G4int n);
  static void SetBarrelStripsPhi(G4int n);
  static void SetAnnularRings(G4int n);
  static void SetAnnularSectors(G4int n);
  static void SetEnableForwardAnnular(G4bool enabled);
  static void SetForwardDistance(G4double d);
  static void SetBackwardDistance(G4double d);
  static void SetDrumModules(G4int n);
  static void SetCapFoldAngle(G4double angle);
  static void SetCapStagger(G4double d);
  static void SetEnableDrum(G4bool enabled);
  static void SetEnableForwardCap(G4bool enabled);
  static void SetEnableBackwardCap(G4bool enabled);

  // ---- Forward-alpha commissioning beam (detector test, default OFF) ----
  // When enabled, PrimaryGeneratorAction replaces the proton beam with alpha
  // particles fired downstream of the target/backing into a forward cone.
  // Detector-response test fixture only; does NOT touch the 11B reaction.
  static G4bool GetForwardAlphaTestBeam();
  static G4double GetForwardAlphaTestBeamEnergy();
  static void SetForwardAlphaTestBeam(G4bool enabled);
  static void SetForwardAlphaTestBeamEnergy(G4double energy);

  // ---- Independent DSSD strip ids from a local hit position ----
  // W1 front face: local-x / azimuthal coordinate, 0..barrelStripsPhi-1.
  static G4int BarrelPhiStripId(G4double local_x, G4double half_x);
  // W1 back face: local-y / polar coordinate, 0..barrelStripsZ-1.
  static G4int BarrelZStripId(G4double local_y, G4double half_y);
  // S3 front face: angular-sector coordinate, 0..annularSectors-1.
  static G4int AnnularSectorStripId(G4double local_x, G4double local_y);
  // S3 back face: radial-ring coordinate, 0..annularRings-1.
  static G4int AnnularRingStripId(G4double local_x, G4double local_y, G4double r_inner, G4double r_outer);

  // Electronic channels per DSSD module.  These are sums of the two faces,
  // not virtual-pixel products.
  static G4int BarrelReadoutChannelCount();
  static G4int AnnularReadoutChannelCount();

  // One-shot startup self-test of the copy_no encoding round-trip.
  static void SelfTestCopyNoEncoding();

  // Non-static command wrappers bound to the /si/ messenger.
  void SetBarrelStripsZCmd(G4int n);
  void SetBarrelStripsPhiCmd(G4int n);
  void SetAnnularRingsCmd(G4int n);
  void SetAnnularSectorsCmd(G4int n);
  void SetEnableForwardAnnularCmd(G4bool enabled);
  void SetForwardDistanceCmd(G4double d);
  void SetBackwardDistanceCmd(G4double d);
  void SetDrumModulesCmd(G4int n);
  void SetCapFoldAngleCmd(G4double angle);
  void SetCapStaggerCmd(G4double d);
  void SetEnableDrumCmd(G4bool enabled);
  void SetEnableForwardCapCmd(G4bool enabled);
  void SetEnableBackwardCapCmd(G4bool enabled);
  void SetForwardAlphaTestBeamCmd(G4bool enabled);
  void SetForwardAlphaTestBeamEnergyCmd(G4double energy);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
