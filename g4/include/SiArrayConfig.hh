#ifndef SiArrayConfig_H
#define SiArrayConfig_H 1

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <memory>

class G4GenericMessenger;

// Centralized, runtime-configurable geometry / strip-segmentation parameters
// for the Si (DSSD-style) array.  Exposed under /si/ so that sensitivity scans
// over strip counts and forward/backward distances can be run without
// recompiling.
//
// IMPORTANT timing note: strip counts are consumed inside SiSD::ProcessHits
// (method (b): a single sensitive module solid + local-coordinate segmentation),
// so they take effect immediately and do NOT require a geometry rebuild.  The
// forward-annular enable flag and the annular axial distances DO affect the
// built geometry; they are read by SiArray::Construct() (executed during
// /run/initialize).  To change them from a macro, set them before
// /run/initialize, or set them and issue /run/reinitializeGeometry.
class SiArrayConfig
{
public:
  SiArrayConfig();
  ~SiArrayConfig();

  // ---- Strip counts (consumed at hit time, no rebuild needed) ----
  static G4int GetBarrelStripsZ();    // strips along the beam axis (local y)
  static G4int GetBarrelStripsPhi();  // strips across the module width (local x)
  static G4int GetAnnularRings();     // radial rings on annular DSSDs
  static G4int GetAnnularSectors();   // angular sectors (spokes) on annular DSSDs

  // ---- Geometry-affecting parameters (consumed at build time) ----
  static G4bool GetEnableForwardAnnular();
  static G4double GetForwardDistance();   // |z - target| of forward annular, Geant4 length units
  static G4double GetBackwardDistance();  // |z - target| of backward annular, Geant4 length units

  static void SetBarrelStripsZ(G4int n);
  static void SetBarrelStripsPhi(G4int n);
  static void SetAnnularRings(G4int n);
  static void SetAnnularSectors(G4int n);
  static void SetEnableForwardAnnular(G4bool enabled);
  static void SetForwardDistance(G4double d);
  static void SetBackwardDistance(G4double d);

  // ---- Forward-alpha commissioning beam (detector test, default OFF) ----
  // When enabled, PrimaryGeneratorAction replaces the proton beam with alpha
  // particles fired downstream of the target/backing into a forward cone that
  // illuminates the forward annular DSSD.  This is a detector-response test
  // fixture only (it does NOT touch the 11B reaction / cross-section / decay
  // model) used to verify that the forward ring registers forward-emitted
  // charged particles; with the default thick target the physics-run forward
  // alphas are absorbed before escaping.
  static G4bool GetForwardAlphaTestBeam();
  static G4double GetForwardAlphaTestBeamEnergy();
  static void SetForwardAlphaTestBeam(G4bool enabled);
  static void SetForwardAlphaTestBeamEnergy(G4double energy);

  // ---- Segment computation from a local hit position ----
  // Barrel module local frame: x across the module width (half_x),
  // y along the beam axis (half_y).  Returns phi_index * n_z + z_index.
  static G4int BarrelSegmentId(G4double local_x, G4double local_y, G4double half_x, G4double half_y);
  // Annular disk local frame: hit lies in the local x-y plane.
  // Returns sector_index * n_rings + ring_index.
  static G4int AnnularSegmentId(G4double local_x, G4double local_y, G4double r_inner, G4double r_outer);

  // Maximum segment id that can be produced for a barrel / annular module with
  // the current configuration (used for range assertions).
  static G4int BarrelSegmentCount();
  static G4int AnnularSegmentCount();

  // One-shot startup self-test: verifies Decode(Encode(...)) round-trips for a
  // representative grid of (type, array, ring, module, segment) tuples and that
  // no field collides.  Aborts via assertion on failure; prints a PASS line.
  static void SelfTestCopyNoEncoding();

  // Non-static command wrappers bound to the /si/ messenger (G4GenericMessenger
  // invokes them on this instance).
  void SetBarrelStripsZCmd(G4int n);
  void SetBarrelStripsPhiCmd(G4int n);
  void SetAnnularRingsCmd(G4int n);
  void SetAnnularSectorsCmd(G4int n);
  void SetEnableForwardAnnularCmd(G4bool enabled);
  void SetForwardDistanceCmd(G4double d);
  void SetBackwardDistanceCmd(G4double d);
  void SetForwardAlphaTestBeamCmd(G4bool enabled);
  void SetForwardAlphaTestBeamEnergyCmd(G4double energy);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
