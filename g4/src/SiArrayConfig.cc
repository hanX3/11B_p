#include "SiArrayConfig.hh"

#include "DetectorChannel.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"

#include <cassert>
#include <cmath>

namespace {
// Defaults follow the all-commercial "connected drum" baseline:
//   * W1 modules (drum + caps): 16 polar x 16 azimuthal strips
//   * S3 annulars: 24 radial rings x 32 angular sectors
G4int g_barrel_strips_z = 16;
G4int g_barrel_strips_phi = 16;
G4int g_annular_rings = 24;
G4int g_annular_sectors = 32;

G4bool g_enable_forward_annular = true;
G4double g_forward_distance = 36. * mm;
G4double g_backward_distance = 36. * mm;

// Connected dodecagonal drum + fish-scale caps.  The drum inscribed radius is
// derived from the module count (see GetDrumInscribedRadius); 12 modules give a
// ~93.74 mm apothem (93.30 mm ideal + face-thickness clearance) and theta
// 75-105.  Caps fold inward 45 deg from the drum end edges
// (theta ~44-75 / 105-136); the outer fish-scale sub-ring sits 8 mm further
// out and half a face period rotated.  8 mm (vs the 4 mm that still let the
// folded inner/outer plates interpenetrate by ~250 um near the beam-axis end)
// clears the overlap checker with ~2 mm of margin over the ~6 mm threshold at
// the default 12 modules / 45 deg fold.
G4int g_drum_modules = 12;
G4double g_cap_fold_angle = 45. * deg;
G4double g_cap_stagger = 8. * mm;
G4bool g_enable_drum = true;
G4bool g_enable_forward_cap = true;
G4bool g_enable_backward_cap = true;

G4bool g_forward_alpha_test_beam = false;
G4double g_forward_alpha_test_beam_energy = 4. * MeV;

// W1 active half-width (50 x 50 mm active area).
constexpr G4double kW1HalfWidth = 25. * mm;
// Half of the 0.5 mm W1 wafer thickness (mirrors SiDetector::map_si_box_par
// "Si_Drum").  The drum faces are flat plates, so the finite thickness drives
// the inner-corner overlap corrected in GetDrumInscribedRadius.
constexpr G4double kW1HalfThickness = 0.25 * mm;
// Perpendicular clearance left between the inner corners of two adjacent drum
// faces and their shared bisector plane, so neighbouring plates do not touch.
constexpr G4double kDrumFaceClearance = 0.05 * mm;

G4int ClampIndex(G4int index, G4int count)
{
  if (index < 0) return 0;
  if (index >= count) return count - 1;
  return index;
}
} // namespace

SiArrayConfig::SiArrayConfig()
{
  DefineCommands();
  SelfTestCopyNoEncoding();
}

SiArrayConfig::~SiArrayConfig() = default;

G4int SiArrayConfig::GetBarrelStripsZ() { return g_barrel_strips_z; }
G4int SiArrayConfig::GetBarrelStripsPhi() { return g_barrel_strips_phi; }
G4int SiArrayConfig::GetAnnularRings() { return g_annular_rings; }
G4int SiArrayConfig::GetAnnularSectors() { return g_annular_sectors; }
G4bool SiArrayConfig::GetEnableForwardAnnular() { return g_enable_forward_annular; }
G4double SiArrayConfig::GetForwardDistance() { return g_forward_distance; }
G4double SiArrayConfig::GetBackwardDistance() { return g_backward_distance; }

G4int SiArrayConfig::GetDrumModules() { return g_drum_modules; }
G4double SiArrayConfig::GetDrumInscribedRadius()
{
  // Regular n-gon of flat, finite-thickness W1 faces.  A zero-thickness face
  // would close exactly at half_width / tan(pi/n) (n=12 -> 93.30 mm), but each
  // plate is kW1HalfThickness*2 thick.  Its *inner* surface (radius apothem -
  // t/2) is a smaller polygon whose faces need less than the full plate width,
  // so the inner corners of neighbouring plates interpenetrate (~57-88 um at
  // n=12).  Push the apothem out by the half-thickness plus a small azimuthal
  // clearance so the inner corners clear the face bisector by kDrumFaceClearance
  // instead of overlapping.  The outer surface opens a matching V-gap that
  // faces away from the target and is harmless.
  const G4double half_angle = CLHEP::pi / static_cast<G4double>(g_drum_modules);
  const G4double apothem = kW1HalfWidth / std::tan(half_angle);
  return apothem + kW1HalfThickness + kDrumFaceClearance / std::sin(half_angle);
}
G4double SiArrayConfig::GetCapFoldAngle() { return g_cap_fold_angle; }
G4double SiArrayConfig::GetCapStagger() { return g_cap_stagger; }
G4bool SiArrayConfig::GetEnableDrum() { return g_enable_drum; }
G4bool SiArrayConfig::GetEnableForwardCap() { return g_enable_forward_cap; }
G4bool SiArrayConfig::GetEnableBackwardCap() { return g_enable_backward_cap; }

void SiArrayConfig::SetBarrelStripsZ(G4int n) { g_barrel_strips_z = n > 0 ? n : 1; }
void SiArrayConfig::SetBarrelStripsPhi(G4int n) { g_barrel_strips_phi = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularRings(G4int n) { g_annular_rings = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularSectors(G4int n) { g_annular_sectors = n > 0 ? n : 1; }
void SiArrayConfig::SetEnableForwardAnnular(G4bool enabled) { g_enable_forward_annular = enabled; }
void SiArrayConfig::SetForwardDistance(G4double d) { g_forward_distance = d > 0. ? d : g_forward_distance; }
void SiArrayConfig::SetBackwardDistance(G4double d) { g_backward_distance = d > 0. ? d : g_backward_distance; }

void SiArrayConfig::SetDrumModules(G4int n)
{
  // Need an even count >= 6: the caps split into two fish-scale sub-rings of
  // n/2, and below 6 faces the inward-folded plates of one sub-ring intersect.
  if (n >= 6 && n % 2 == 0) g_drum_modules = n;
}
void SiArrayConfig::SetCapFoldAngle(G4double angle)
{
  if (angle > 0. && angle < 90. * deg) g_cap_fold_angle = angle;
}
void SiArrayConfig::SetCapStagger(G4double d) { g_cap_stagger = d > 0. ? d : g_cap_stagger; }
void SiArrayConfig::SetEnableDrum(G4bool enabled) { g_enable_drum = enabled; }
void SiArrayConfig::SetEnableForwardCap(G4bool enabled) { g_enable_forward_cap = enabled; }
void SiArrayConfig::SetEnableBackwardCap(G4bool enabled) { g_enable_backward_cap = enabled; }

G4bool SiArrayConfig::GetForwardAlphaTestBeam() { return g_forward_alpha_test_beam; }
G4double SiArrayConfig::GetForwardAlphaTestBeamEnergy() { return g_forward_alpha_test_beam_energy; }
void SiArrayConfig::SetForwardAlphaTestBeam(G4bool enabled) { g_forward_alpha_test_beam = enabled; }
void SiArrayConfig::SetForwardAlphaTestBeamEnergy(G4double energy) { g_forward_alpha_test_beam_energy = energy > 0. ? energy : g_forward_alpha_test_beam_energy; }

G4int SiArrayConfig::BarrelReadoutChannelCount() { return g_barrel_strips_phi + g_barrel_strips_z; }
G4int SiArrayConfig::AnnularReadoutChannelCount() { return g_annular_sectors + g_annular_rings; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::BarrelPhiStripId(G4double local_x, G4double half_x)
{
  const G4double span_x = 2. * half_x;
  return span_x > 0. ? ClampIndex(static_cast<G4int>((local_x + half_x) / span_x * g_barrel_strips_phi),
                                    g_barrel_strips_phi)
                     : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::BarrelZStripId(G4double local_y, G4double half_y)
{
  const G4double span_y = 2. * half_y;
  return span_y > 0. ? ClampIndex(static_cast<G4int>((local_y + half_y) / span_y * g_barrel_strips_z),
                                    g_barrel_strips_z)
                     : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::AnnularSectorStripId(G4double local_x, G4double local_y)
{
  G4double phi = std::atan2(local_y, local_x); // [-pi, pi]
  if (phi < 0.) phi += CLHEP::twopi;
  return ClampIndex(static_cast<G4int>(phi / CLHEP::twopi * g_annular_sectors), g_annular_sectors);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::AnnularRingStripId(G4double local_x, G4double local_y, G4double r_inner, G4double r_outer)
{
  const G4double r = std::hypot(local_x, local_y);
  const G4double span_r = r_outer - r_inner;
  return span_r > 0. ? ClampIndex(static_cast<G4int>((r - r_inner) / span_r * g_annular_rings),
                                    g_annular_rings)
                     : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::SelfTestCopyNoEncoding()
{
  const G4int types[] = {1, 2, 3};
  const G4int arrays[] = {0};
  const G4int rings[] = {1, 2, 3, 4, 5, 42, kMaxRingId};
  const G4int modules[] = {0, 1, 7, 11, 15, kMaxModuleId};
  const G4int segments[] = {0, 1, 9, 15, 255, 767, kMaxSegmentId};

  G4int checked = 0;
  for (G4int type : types) {
    for (G4int array : arrays) {
      for (G4int ring : rings) {
        for (G4int module : modules) {
          for (G4int segment : segments) {
            const G4int copy_no = EncodeDetectorCopyNo(static_cast<DetectorType>(type), array, ring, module, segment);
            const DetectorChannel back = DecodeDetectorCopyNo(copy_no);
            assert(back.detector_type == type && "copy_no type round-trip failed");
            assert(back.array_id == array && "copy_no array round-trip failed");
            assert(back.ring_id == ring && "copy_no ring round-trip failed");
            assert(back.module_id == module && "copy_no module round-trip failed");
            assert(back.segment_id == segment && "copy_no segment round-trip failed");
            assert(copy_no / kTypeUnit == type && "copy_no type not top field");
            ++checked;
          }
        }
      }
    }
  }
  G4cout << "[SiArrayConfig] copy_no encoding self-test PASS (" << checked << " tuples reversible)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::SetBarrelStripsZCmd(G4int n) { SetBarrelStripsZ(n); }
void SiArrayConfig::SetBarrelStripsPhiCmd(G4int n) { SetBarrelStripsPhi(n); }
void SiArrayConfig::SetAnnularRingsCmd(G4int n) { SetAnnularRings(n); }
void SiArrayConfig::SetAnnularSectorsCmd(G4int n) { SetAnnularSectors(n); }
void SiArrayConfig::SetEnableForwardAnnularCmd(G4bool enabled) { SetEnableForwardAnnular(enabled); }
void SiArrayConfig::SetForwardDistanceCmd(G4double d) { SetForwardDistance(d); }
void SiArrayConfig::SetBackwardDistanceCmd(G4double d) { SetBackwardDistance(d); }
void SiArrayConfig::SetDrumModulesCmd(G4int n) { SetDrumModules(n); }
void SiArrayConfig::SetCapFoldAngleCmd(G4double angle) { SetCapFoldAngle(angle); }
void SiArrayConfig::SetCapStaggerCmd(G4double d) { SetCapStagger(d); }
void SiArrayConfig::SetEnableDrumCmd(G4bool enabled) { SetEnableDrum(enabled); }
void SiArrayConfig::SetEnableForwardCapCmd(G4bool enabled) { SetEnableForwardCap(enabled); }
void SiArrayConfig::SetEnableBackwardCapCmd(G4bool enabled) { SetEnableBackwardCap(enabled); }
void SiArrayConfig::SetForwardAlphaTestBeamCmd(G4bool enabled) { SetForwardAlphaTestBeam(enabled); }
void SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd(G4double energy) { SetForwardAlphaTestBeamEnergy(energy); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::PrintConfigCommand()
{
  G4cout << "\n===== Si array runtime configuration (connected drum baseline) =====" << G4endl
         << "  barrelStripsZ (W1 polar strips)   = " << g_barrel_strips_z << G4endl
         << "  barrelStripsPhi (W1 azim. strips) = " << g_barrel_strips_phi << G4endl
         << "  annularRings (S3)                 = " << g_annular_rings << G4endl
         << "  annularSectors (S3)               = " << g_annular_sectors << G4endl
         << "  enableForwardAnnular              = " << (g_enable_forward_annular ? "true" : "false") << G4endl
         << "  forwardDistance                   = " << g_forward_distance / mm << " mm" << G4endl
         << "  backwardDistance                  = " << g_backward_distance / mm << " mm" << G4endl
         << "  drumModules                       = " << g_drum_modules << G4endl
         << "  drum inscribed radius (derived)   = " << GetDrumInscribedRadius() / mm << " mm" << G4endl
         << "  capFoldAngle                      = " << g_cap_fold_angle / deg << " deg" << G4endl
         << "  capStagger                        = " << g_cap_stagger / mm << " mm" << G4endl
         << "  enableDrum                        = " << (g_enable_drum ? "true" : "false") << G4endl
         << "  enableForwardCap                  = " << (g_enable_forward_cap ? "true" : "false") << G4endl
         << "  enableBackwardCap                 = " << (g_enable_backward_cap ? "true" : "false") << G4endl
         << "  forwardAlphaTestBeam              = " << (g_forward_alpha_test_beam ? "true" : "false")
         << " (" << g_forward_alpha_test_beam_energy / MeV << " MeV)" << G4endl
         << "  W1 readout channels/mod           = " << BarrelReadoutChannelCount()
         << " (" << g_barrel_strips_phi << " front + " << g_barrel_strips_z << " back)" << G4endl
         << "  S3 readout channels/mod           = " << AnnularReadoutChannelCount()
         << " (" << g_annular_sectors << " front + " << g_annular_rings << " back)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/si/", "Si array geometry / strip-segmentation controls");

  messenger->DeclareMethod("barrelStripsZ", &SiArrayConfig::SetBarrelStripsZCmd,
                           "W1 strips, polar direction (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("barrelStripsPhi", &SiArrayConfig::SetBarrelStripsPhiCmd,
                           "W1 strips, azimuthal direction (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("annularRings", &SiArrayConfig::SetAnnularRingsCmd,
                           "S3 annular radial rings (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("annularSectors", &SiArrayConfig::SetAnnularSectorsCmd,
                           "S3 annular angular sectors (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("enableForwardAnnular", &SiArrayConfig::SetEnableForwardAnnularCmd,
                           "Enable the forward S3 annular (geometry; set before /run/initialize)");
  messenger->DeclareMethodWithUnit("forwardDistance", "mm", &SiArrayConfig::SetForwardDistanceCmd,
                                   "Forward annular |z-target| distance (geometry)");
  messenger->DeclareMethodWithUnit("backwardDistance", "mm", &SiArrayConfig::SetBackwardDistanceCmd,
                                   "Backward annular |z-target| distance (geometry)");
  messenger->DeclareMethod("drumModules", &SiArrayConfig::SetDrumModulesCmd,
                           "W1 faces of the drum, even and >= 6; inscribed radius follows (geometry)");
  messenger->DeclareMethodWithUnit("capFoldAngle", "deg", &SiArrayConfig::SetCapFoldAngleCmd,
                                   "Inward fold of the cap plates from the drum plane, (0,90) deg (geometry)");
  messenger->DeclareMethodWithUnit("capStagger", "mm", &SiArrayConfig::SetCapStaggerCmd,
                                   "Radial offset of the outer fish-scale cap sub-ring (geometry)");
  messenger->DeclareMethod("enableDrum", &SiArrayConfig::SetEnableDrumCmd,
                           "Enable the equatorial W1 drum (geometry)");
  messenger->DeclareMethod("enableForwardCap", &SiArrayConfig::SetEnableForwardCapCmd,
                           "Enable the forward W1 fish-scale cap (geometry)");
  messenger->DeclareMethod("enableBackwardCap", &SiArrayConfig::SetEnableBackwardCapCmd,
                           "Enable the backward W1 fish-scale cap (geometry)");
  messenger->DeclareMethod("forwardAlphaTestBeam", &SiArrayConfig::SetForwardAlphaTestBeamCmd,
                           "Detector test only: fire a forward alpha beam (default off)");
  messenger->DeclareMethodWithUnit("forwardAlphaTestBeamEnergy", "MeV", &SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd,
                                   "Kinetic energy of the forward-alpha commissioning beam");
  messenger->DeclareMethod("printConfig", &SiArrayConfig::PrintConfigCommand,
                           "Print current Si array configuration");
}
