#include "SiArrayConfig.hh"

#include "DetectorChannel.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <cassert>
#include <cmath>

namespace {
// Defaults follow the task specification: barrel 16 strips along the beam,
// single-sided (phi=1); annular DSSDs 16 radial rings x 16 angular sectors.
G4int g_barrel_strips_z = 16;
G4int g_barrel_strips_phi = 1;
G4int g_annular_rings = 16;
G4int g_annular_sectors = 16;

G4bool g_enable_forward_annular = true;
G4double g_forward_distance = 120. * mm;
G4double g_backward_distance = 120. * mm;

G4bool g_forward_alpha_test_beam = false;
G4double g_forward_alpha_test_beam_energy = 4. * MeV;

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

void SiArrayConfig::SetBarrelStripsZ(G4int n) { g_barrel_strips_z = n > 0 ? n : 1; }
void SiArrayConfig::SetBarrelStripsPhi(G4int n) { g_barrel_strips_phi = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularRings(G4int n) { g_annular_rings = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularSectors(G4int n) { g_annular_sectors = n > 0 ? n : 1; }
void SiArrayConfig::SetEnableForwardAnnular(G4bool enabled) { g_enable_forward_annular = enabled; }
void SiArrayConfig::SetForwardDistance(G4double d) { g_forward_distance = d > 0. ? d : g_forward_distance; }
void SiArrayConfig::SetBackwardDistance(G4double d) { g_backward_distance = d > 0. ? d : g_backward_distance; }

G4bool SiArrayConfig::GetForwardAlphaTestBeam() { return g_forward_alpha_test_beam; }
G4double SiArrayConfig::GetForwardAlphaTestBeamEnergy() { return g_forward_alpha_test_beam_energy; }
void SiArrayConfig::SetForwardAlphaTestBeam(G4bool enabled) { g_forward_alpha_test_beam = enabled; }
void SiArrayConfig::SetForwardAlphaTestBeamEnergy(G4double energy) { g_forward_alpha_test_beam_energy = energy > 0. ? energy : g_forward_alpha_test_beam_energy; }

G4int SiArrayConfig::BarrelSegmentCount() { return g_barrel_strips_phi * g_barrel_strips_z; }
G4int SiArrayConfig::AnnularSegmentCount() { return g_annular_sectors * g_annular_rings; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::BarrelSegmentId(G4double local_x, G4double local_y, G4double half_x, G4double half_y)
{
  const G4int n_z = g_barrel_strips_z;
  const G4int n_phi = g_barrel_strips_phi;

  const G4double span_x = 2. * half_x;
  const G4double span_y = 2. * half_y;
  const G4int z_index = span_y > 0. ? ClampIndex(static_cast<G4int>((local_y + half_y) / span_y * n_z), n_z) : 0;
  const G4int phi_index = span_x > 0. ? ClampIndex(static_cast<G4int>((local_x + half_x) / span_x * n_phi), n_phi) : 0;

  return phi_index * n_z + z_index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int SiArrayConfig::AnnularSegmentId(G4double local_x, G4double local_y, G4double r_inner, G4double r_outer)
{
  const G4int n_rings = g_annular_rings;
  const G4int n_sectors = g_annular_sectors;

  const G4double r = std::hypot(local_x, local_y);
  const G4double span_r = r_outer - r_inner;
  const G4int ring_index = span_r > 0. ? ClampIndex(static_cast<G4int>((r - r_inner) / span_r * n_rings), n_rings) : 0;

  G4double phi = std::atan2(local_y, local_x); // [-pi, pi]
  if (phi < 0.) phi += CLHEP::twopi;
  const G4int sector_index = ClampIndex(static_cast<G4int>(phi / CLHEP::twopi * n_sectors), n_sectors);

  return sector_index * n_rings + ring_index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::SelfTestCopyNoEncoding()
{
  const G4int types[] = {1, 2, 3};
  const G4int arrays[] = {0};
  const G4int rings[] = {1, 2, 3, 42, kMaxRingId};
  const G4int modules[] = {0, 1, 7, 15, kMaxModuleId};
  const G4int segments[] = {0, 1, 9, 15, 255, kMaxSegmentId};

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
            // detector_type must remain recoverable as copy_no / kTypeUnit, as
            // relied on by the validation scripts.
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
void SiArrayConfig::SetForwardAlphaTestBeamCmd(G4bool enabled) { SetForwardAlphaTestBeam(enabled); }
void SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd(G4double energy) { SetForwardAlphaTestBeamEnergy(energy); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::PrintConfigCommand()
{
  G4cout << "\n===== Si array runtime configuration =====" << G4endl
         << "  barrelStripsZ        = " << g_barrel_strips_z << G4endl
         << "  barrelStripsPhi      = " << g_barrel_strips_phi << G4endl
         << "  annularRings         = " << g_annular_rings << G4endl
         << "  annularSectors       = " << g_annular_sectors << G4endl
         << "  enableForwardAnnular = " << (g_enable_forward_annular ? "true" : "false") << G4endl
         << "  forwardDistance      = " << g_forward_distance / mm << " mm" << G4endl
         << "  backwardDistance     = " << g_backward_distance / mm << " mm" << G4endl
         << "  forwardAlphaTestBeam = " << (g_forward_alpha_test_beam ? "true" : "false")
         << " (" << g_forward_alpha_test_beam_energy / MeV << " MeV)" << G4endl
         << "  barrel segments/mod  = " << BarrelSegmentCount() << G4endl
         << "  annular segments/mod = " << AnnularSegmentCount() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArrayConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/si/", "Si array geometry / strip-segmentation controls");

  messenger->DeclareMethod("barrelStripsZ", &SiArrayConfig::SetBarrelStripsZCmd,
                           "Barrel strips along the beam axis (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("barrelStripsPhi", &SiArrayConfig::SetBarrelStripsPhiCmd,
                           "Barrel strips across the module width (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("annularRings", &SiArrayConfig::SetAnnularRingsCmd,
                           "Annular DSSD radial rings (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("annularSectors", &SiArrayConfig::SetAnnularSectorsCmd,
                           "Annular DSSD angular sectors (takes effect immediately, no rebuild)");
  messenger->DeclareMethod("enableForwardAnnular", &SiArrayConfig::SetEnableForwardAnnularCmd,
                           "Enable the forward annular Si (geometry; set before /run/initialize)");
  messenger->DeclareMethodWithUnit("forwardDistance", "mm", &SiArrayConfig::SetForwardDistanceCmd,
                                   "Forward annular |z-target| distance (geometry; set before /run/initialize)");
  messenger->DeclareMethodWithUnit("backwardDistance", "mm", &SiArrayConfig::SetBackwardDistanceCmd,
                                   "Backward annular |z-target| distance (geometry; set before /run/initialize)");
  messenger->DeclareMethod("forwardAlphaTestBeam", &SiArrayConfig::SetForwardAlphaTestBeamCmd,
                           "Detector test only: fire a forward alpha beam onto the forward annular (default off)");
  messenger->DeclareMethodWithUnit("forwardAlphaTestBeamEnergy", "MeV", &SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd,
                                   "Kinetic energy of the forward-alpha commissioning beam");
  messenger->DeclareMethod("printConfig", &SiArrayConfig::PrintConfigCommand,
                           "Print current Si array configuration");
}
