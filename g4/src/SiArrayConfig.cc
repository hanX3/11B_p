#include "SiArrayConfig.hh"

#include "DetectorChannel.hh"

#include "G4GenericMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace {
// W1 drum: real W1 16 x 16 electronic strips.
G4int g_barrel_strips_z = 16;
G4int g_barrel_strips_phi = 16;

// S3: 24 radial rings x 32 sectors.
G4int g_annular_rings = 24;
G4int g_annular_sectors = 32;

// CAKE/MMM-inspired wedge: 16 radial rings x 8 sectors.
G4int g_lampshade_rings = 16;
G4int g_lampshade_sectors = 8;

G4bool g_enable_forward_annular = true;
G4bool g_enable_backward_annular = true;

// The S3 outer radius (35 mm) lies on the same ray as the lampshade narrow
// edge.  With the default lampshade geometry this is 159.291607 mm.
G4double g_forward_distance = 159.2916073617 * mm;
G4double g_backward_distance = 159.2916073617 * mm;

G4int g_drum_modules = 12;
G4bool g_enable_drum = true;

// Twelvefold CAKE-like lampshades.  Each wedge centre is 98 mm from the target
// at theta=40 degrees (forward) or 140 degrees (backward).  The 102.5-mm long
// detector is tangent to the target-centred sphere at its centre, giving
// narrow/wide edge angles of approximately 12.392/67.608 degrees.  Module
// centres are spaced by 30 degrees and share the same phi centres as the
// 12-module W1 drum.  Each wedge spans 28 degrees, leaving a regular 2-degree
// inter-module gap.
G4bool g_enable_forward_cap = true;
G4bool g_enable_backward_cap = true;
G4int g_lampshade_modules = 12;
G4double g_lampshade_center_distance = 98. * mm;
G4double g_lampshade_center_angle = 40. * deg;
G4double g_lampshade_length = 102.5 * mm;
G4double g_lampshade_module_span = 28. * deg;
G4double g_lampshade_thickness = 0.4 * mm;

G4bool g_forward_alpha_test_beam = false;
G4double g_forward_alpha_test_beam_energy = 4. * MeV;

G4bool g_fixed_alpha_test_beam = false;
G4double g_fixed_alpha_test_beam_energy = 5. * MeV;
G4double g_fixed_alpha_test_beam_theta = 171.607793103 * deg;
G4double g_fixed_alpha_test_beam_phi = 5.625 * deg;

constexpr G4double kW1HalfWidth = 25. * mm;
constexpr G4double kW1HalfThickness = 0.25 * mm;
constexpr G4double kDrumFaceClearance = 0.05 * mm;

G4int ClampIndex(G4int index, G4int count)
{
  if (count <= 1) return 0;
  if (index < 0) return 0;
  if (index >= count) return count - 1;
  return index;
}

struct LampshadeEdges
{
  G4double inner_radius = 0.;
  G4double inner_distance = 0.;
  G4double outer_radius = 0.;
  G4double outer_distance = 0.;
};

LampshadeEdges CalculateLampshadeEdges(G4double centre_distance,
                                       G4double centre_angle,
                                       G4double length)
{
  const G4double half_length = 0.5 * length;
  const G4double sin_theta = std::sin(centre_angle);
  const G4double cos_theta = std::cos(centre_angle);

  // The detector centre is D*e_r and its long axis follows +e_theta.  Hence
  // narrow edge = centre - L/2*e_theta and wide edge = centre + L/2*e_theta.
  return {
      centre_distance * sin_theta - half_length * cos_theta,
      centre_distance * cos_theta + half_length * sin_theta,
      centre_distance * sin_theta + half_length * cos_theta,
      centre_distance * cos_theta - half_length * sin_theta};
}

LampshadeEdges CalculateLampshadeEdges()
{
  return CalculateLampshadeEdges(g_lampshade_center_distance,
                                 g_lampshade_center_angle,
                                 g_lampshade_length);
}

G4bool LampshadeGeometryIsValid(G4double centre_distance,
                                G4double centre_angle,
                                G4double length)
{
  if (centre_distance <= 0. || length <= 0.) return false;
  if (centre_angle <= 0. || centre_angle >= 90. * deg) return false;

  const LampshadeEdges edges = CalculateLampshadeEdges(centre_distance,
                                                       centre_angle,
                                                       length);
  return edges.inner_radius > 0.
         && edges.inner_distance > 0.
         && edges.outer_radius > edges.inner_radius
         && edges.outer_distance > 0.;
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
G4int SiArrayConfig::GetLampshadeRings() { return g_lampshade_rings; }
G4int SiArrayConfig::GetLampshadeSectors() { return g_lampshade_sectors; }

G4bool SiArrayConfig::GetEnableForwardAnnular() { return g_enable_forward_annular; }
G4bool SiArrayConfig::GetEnableBackwardAnnular() { return g_enable_backward_annular; }
G4double SiArrayConfig::GetForwardDistance() { return g_forward_distance; }
G4double SiArrayConfig::GetBackwardDistance() { return g_backward_distance; }

G4int SiArrayConfig::GetDrumModules() { return g_drum_modules; }
G4double SiArrayConfig::GetDrumInscribedRadius()
{
  const G4double half_angle = CLHEP::pi / static_cast<G4double>(g_drum_modules);
  const G4double apothem = kW1HalfWidth / std::tan(half_angle);
  return apothem + kW1HalfThickness
         + kDrumFaceClearance / std::sin(half_angle);
}
G4bool SiArrayConfig::GetEnableDrum() { return g_enable_drum; }

G4bool SiArrayConfig::GetEnableForwardCap() { return g_enable_forward_cap; }
G4bool SiArrayConfig::GetEnableBackwardCap() { return g_enable_backward_cap; }
G4int SiArrayConfig::GetLampshadeModules() { return g_lampshade_modules; }
G4double SiArrayConfig::GetLampshadeCenterDistance() { return g_lampshade_center_distance; }
G4double SiArrayConfig::GetLampshadeCenterAngle() { return g_lampshade_center_angle; }
G4double SiArrayConfig::GetLampshadeLength() { return g_lampshade_length; }
G4double SiArrayConfig::GetLampshadeModuleSpan() { return g_lampshade_module_span; }
G4double SiArrayConfig::GetLampshadeThickness() { return g_lampshade_thickness; }
G4double SiArrayConfig::GetLampshadeInnerRadius() { return CalculateLampshadeEdges().inner_radius; }
G4double SiArrayConfig::GetLampshadeInnerDistance() { return CalculateLampshadeEdges().inner_distance; }
G4double SiArrayConfig::GetLampshadeOuterRadius() { return CalculateLampshadeEdges().outer_radius; }
G4double SiArrayConfig::GetLampshadeOuterDistance() { return CalculateLampshadeEdges().outer_distance; }
G4double SiArrayConfig::GetLampshadeInnerAngle()
{
  return std::atan2(GetLampshadeInnerRadius(), GetLampshadeInnerDistance());
}
G4double SiArrayConfig::GetLampshadeOuterAngle()
{
  return std::atan2(GetLampshadeOuterRadius(), GetLampshadeOuterDistance());
}
G4double SiArrayConfig::GetLampshadeInnerHalfWidth()
{
  return GetLampshadeInnerRadius() * std::tan(0.5 * g_lampshade_module_span);
}
G4double SiArrayConfig::GetLampshadeOuterHalfWidth()
{
  return GetLampshadeOuterRadius() * std::tan(0.5 * g_lampshade_module_span);
}
G4double SiArrayConfig::GetLampshadeHoleApothem() { return GetLampshadeInnerRadius(); }

// Obsolete W1-cap controls retained only for old macro compatibility.
G4double SiArrayConfig::GetCapFoldAngle() { return 0.; }
G4double SiArrayConfig::GetCapClearance() { return 0.; }
G4double SiArrayConfig::GetCapStagger() { return 0.; }

void SiArrayConfig::SetBarrelStripsZ(G4int n) { g_barrel_strips_z = n > 0 ? n : 1; }
void SiArrayConfig::SetBarrelStripsPhi(G4int n) { g_barrel_strips_phi = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularRings(G4int n) { g_annular_rings = n > 0 ? n : 1; }
void SiArrayConfig::SetAnnularSectors(G4int n) { g_annular_sectors = n > 0 ? n : 1; }
void SiArrayConfig::SetLampshadeRings(G4int n) { g_lampshade_rings = n > 0 ? n : 1; }
void SiArrayConfig::SetLampshadeSectors(G4int n) { g_lampshade_sectors = n > 0 ? n : 1; }

void SiArrayConfig::SetEnableForwardAnnular(G4bool enabled) { g_enable_forward_annular = enabled; }
void SiArrayConfig::SetEnableBackwardAnnular(G4bool enabled) { g_enable_backward_annular = enabled; }
void SiArrayConfig::SetForwardDistance(G4double d)
{
  if (d > 0.) g_forward_distance = d;
}
void SiArrayConfig::SetBackwardDistance(G4double d)
{
  if (d > 0.) g_backward_distance = d;
}

void SiArrayConfig::SetDrumModules(G4int n)
{
  if (n >= 6) g_drum_modules = n;
}
void SiArrayConfig::SetEnableDrum(G4bool enabled) { g_enable_drum = enabled; }

void SiArrayConfig::SetEnableForwardCap(G4bool enabled) { g_enable_forward_cap = enabled; }
void SiArrayConfig::SetEnableBackwardCap(G4bool enabled) { g_enable_backward_cap = enabled; }
void SiArrayConfig::SetLampshadeModules(G4int n)
{
  if (n < 3) return;
  const G4double cell_span = CLHEP::twopi / static_cast<G4double>(n);
  if (g_lampshade_module_span < cell_span) g_lampshade_modules = n;
}
void SiArrayConfig::SetLampshadeCenterDistance(G4double d)
{
  if (LampshadeGeometryIsValid(d,
                               g_lampshade_center_angle,
                               g_lampshade_length)) {
    g_lampshade_center_distance = d;
  }
}
void SiArrayConfig::SetLampshadeCenterAngle(G4double angle)
{
  if (LampshadeGeometryIsValid(g_lampshade_center_distance,
                               angle,
                               g_lampshade_length)) {
    g_lampshade_center_angle = angle;
  }
}
void SiArrayConfig::SetLampshadeLength(G4double d)
{
  if (LampshadeGeometryIsValid(g_lampshade_center_distance,
                               g_lampshade_center_angle,
                               d)) {
    g_lampshade_length = d;
  }
}
void SiArrayConfig::SetLampshadeModuleSpan(G4double angle)
{
  const G4double cell_span = CLHEP::twopi / static_cast<G4double>(g_lampshade_modules);
  if (angle > 0. && angle < cell_span) g_lampshade_module_span = angle;
}
void SiArrayConfig::SetLampshadeThickness(G4double d)
{
  if (d > 0.) g_lampshade_thickness = d;
}

void SiArrayConfig::SetLampshadeHoleApothem(G4double)
{
  G4cout << "[SiArrayConfig] /si/lampshadeHoleApothem is obsolete: use "
         << "/si/lampshadeCenterDistance and /si/lampshadeCenterAngle." << G4endl;
}
void SiArrayConfig::SetLampshadeInnerDistance(G4double)
{
  G4cout << "[SiArrayConfig] /si/lampshadeInnerDistance is obsolete: use "
         << "/si/lampshadeCenterDistance and /si/lampshadeCenterAngle." << G4endl;
}
void SiArrayConfig::SetLampshadeOuterAngle(G4double)
{
  G4cout << "[SiArrayConfig] /si/lampshadeOuterAngle is obsolete: use "
         << "/si/lampshadeCenterAngle." << G4endl;
}

void SiArrayConfig::SetCapFoldAngle(G4double)
{
  G4cout << "[SiArrayConfig] /si/capFoldAngle is obsolete: square W1 caps were replaced by wedge lampshades." << G4endl;
}
void SiArrayConfig::SetCapClearance(G4double)
{
  G4cout << "[SiArrayConfig] /si/capClearance is obsolete: use /si/lampshadeModuleSpan." << G4endl;
}
void SiArrayConfig::SetCapStagger(G4double)
{
  G4cout << "[SiArrayConfig] /si/capStagger is obsolete: no fish-scale W1 caps remain." << G4endl;
}

G4bool SiArrayConfig::GetForwardAlphaTestBeam() { return g_forward_alpha_test_beam; }
G4double SiArrayConfig::GetForwardAlphaTestBeamEnergy() { return g_forward_alpha_test_beam_energy; }
void SiArrayConfig::SetForwardAlphaTestBeam(G4bool enabled) { g_forward_alpha_test_beam = enabled; }
void SiArrayConfig::SetForwardAlphaTestBeamEnergy(G4double energy)
{
  if (energy > 0.) g_forward_alpha_test_beam_energy = energy;
}

G4bool SiArrayConfig::GetFixedAlphaTestBeam() { return g_fixed_alpha_test_beam; }
G4double SiArrayConfig::GetFixedAlphaTestBeamEnergy() { return g_fixed_alpha_test_beam_energy; }
G4double SiArrayConfig::GetFixedAlphaTestBeamTheta() { return g_fixed_alpha_test_beam_theta; }
G4double SiArrayConfig::GetFixedAlphaTestBeamPhi() { return g_fixed_alpha_test_beam_phi; }
void SiArrayConfig::SetFixedAlphaTestBeam(G4bool enabled) { g_fixed_alpha_test_beam = enabled; }
void SiArrayConfig::SetFixedAlphaTestBeamEnergy(G4double energy)
{
  if (energy > 0.) g_fixed_alpha_test_beam_energy = energy;
}
void SiArrayConfig::SetFixedAlphaTestBeamTheta(G4double theta)
{
  if (theta >= 0. && theta <= CLHEP::pi) g_fixed_alpha_test_beam_theta = theta;
}
void SiArrayConfig::SetFixedAlphaTestBeamPhi(G4double phi)
{
  g_fixed_alpha_test_beam_phi = std::fmod(phi, CLHEP::twopi);
  if (g_fixed_alpha_test_beam_phi < 0.) g_fixed_alpha_test_beam_phi += CLHEP::twopi;
}

G4int SiArrayConfig::BarrelPhiStripId(G4double local_x, G4double half_x)
{
  const G4double span_x = 2. * half_x;
  return span_x > 0.
             ? ClampIndex(static_cast<G4int>((local_x + half_x) / span_x * g_barrel_strips_phi),
                          g_barrel_strips_phi)
             : 0;
}

G4int SiArrayConfig::BarrelZStripId(G4double local_y, G4double half_y)
{
  const G4double span_y = 2. * half_y;
  return span_y > 0.
             ? ClampIndex(static_cast<G4int>((local_y + half_y) / span_y * g_barrel_strips_z),
                          g_barrel_strips_z)
             : 0;
}

G4int SiArrayConfig::AnnularSectorStripId(G4double local_x, G4double local_y)
{
  G4double phi = std::atan2(local_y, local_x);
  if (phi < 0.) phi += CLHEP::twopi;
  return ClampIndex(static_cast<G4int>(phi / CLHEP::twopi * g_annular_sectors),
                    g_annular_sectors);
}

G4int SiArrayConfig::AnnularRingStripId(G4double local_x,
                                        G4double local_y,
                                        G4double r_inner,
                                        G4double r_outer)
{
  const G4double radius = std::hypot(local_x, local_y);
  const G4double span = r_outer - r_inner;
  return span > 0.
             ? ClampIndex(static_cast<G4int>((radius - r_inner) / span * g_annular_rings),
                          g_annular_rings)
             : 0;
}

G4double SiArrayConfig::LampshadeHalfWidthAtY(G4double local_y,
                                              G4double half_length,
                                              G4double inner_half_width,
                                              G4double outer_half_width)
{
  if (half_length <= 0.) return inner_half_width;
  const G4double fraction = std::clamp((local_y + half_length) / (2. * half_length), 0., 1.);
  return inner_half_width + fraction * (outer_half_width - inner_half_width);
}

G4int SiArrayConfig::LampshadeSectorStripId(G4double local_x,
                                            G4double local_y,
                                            G4double half_length,
                                            G4double inner_half_width,
                                            G4double outer_half_width)
{
  const G4double half_width = LampshadeHalfWidthAtY(local_y,
                                                    half_length,
                                                    inner_half_width,
                                                    outer_half_width);
  return half_width > 0.
             ? ClampIndex(static_cast<G4int>((local_x + half_width) / (2. * half_width)
                                             * g_lampshade_sectors),
                          g_lampshade_sectors)
             : 0;
}

G4int SiArrayConfig::LampshadeRingStripId(G4double local_y, G4double half_length)
{
  return half_length > 0.
             ? ClampIndex(static_cast<G4int>((local_y + half_length) / (2. * half_length)
                                             * g_lampshade_rings),
                          g_lampshade_rings)
             : 0;
}

G4int SiArrayConfig::BarrelReadoutChannelCount()
{
  return g_barrel_strips_phi + g_barrel_strips_z;
}
G4int SiArrayConfig::AnnularReadoutChannelCount()
{
  return g_annular_sectors + g_annular_rings;
}
G4int SiArrayConfig::LampshadeReadoutChannelCount()
{
  return g_lampshade_sectors + g_lampshade_rings;
}

void SiArrayConfig::SelfTestCopyNoEncoding()
{
  const G4int types[] = {1, 2, 3};
  const G4int arrays[] = {0};
  const G4int rings[] = {1, 2, 3, 4, 5, 42, kMaxRingId};
  const G4int modules[] = {0, 1, 5, 7, 11, kMaxModuleId};
  const G4int segments[] = {0, 1, 7, 15, 31, 255, kMaxSegmentId};

  G4int checked = 0;
  for (G4int type : types) {
    for (G4int array : arrays) {
      for (G4int ring : rings) {
        for (G4int module : modules) {
          for (G4int segment : segments) {
            const G4int copy_no = EncodeDetectorCopyNo(static_cast<DetectorType>(type),
                                                       array,
                                                       ring,
                                                       module,
                                                       segment);
            const DetectorChannel back = DecodeDetectorCopyNo(copy_no);
            assert(back.detector_type == type);
            assert(back.array_id == array);
            assert(back.ring_id == ring);
            assert(back.module_id == module);
            assert(back.segment_id == segment);
            ++checked;
          }
        }
      }
    }
  }
  G4cout << "[SiArrayConfig] copy_no encoding self-test PASS ("
         << checked << " tuples reversible)" << G4endl;
}

void SiArrayConfig::SetBarrelStripsZCmd(G4int n) { SetBarrelStripsZ(n); }
void SiArrayConfig::SetBarrelStripsPhiCmd(G4int n) { SetBarrelStripsPhi(n); }
void SiArrayConfig::SetAnnularRingsCmd(G4int n) { SetAnnularRings(n); }
void SiArrayConfig::SetAnnularSectorsCmd(G4int n) { SetAnnularSectors(n); }
void SiArrayConfig::SetLampshadeRingsCmd(G4int n) { SetLampshadeRings(n); }
void SiArrayConfig::SetLampshadeSectorsCmd(G4int n) { SetLampshadeSectors(n); }
void SiArrayConfig::SetEnableForwardAnnularCmd(G4bool enabled) { SetEnableForwardAnnular(enabled); }
void SiArrayConfig::SetEnableBackwardAnnularCmd(G4bool enabled) { SetEnableBackwardAnnular(enabled); }
void SiArrayConfig::SetForwardDistanceCmd(G4double d) { SetForwardDistance(d); }
void SiArrayConfig::SetBackwardDistanceCmd(G4double d) { SetBackwardDistance(d); }
void SiArrayConfig::SetDrumModulesCmd(G4int n) { SetDrumModules(n); }
void SiArrayConfig::SetEnableDrumCmd(G4bool enabled) { SetEnableDrum(enabled); }
void SiArrayConfig::SetEnableForwardCapCmd(G4bool enabled) { SetEnableForwardCap(enabled); }
void SiArrayConfig::SetEnableBackwardCapCmd(G4bool enabled) { SetEnableBackwardCap(enabled); }
void SiArrayConfig::SetLampshadeModulesCmd(G4int n) { SetLampshadeModules(n); }
void SiArrayConfig::SetLampshadeCenterDistanceCmd(G4double d) { SetLampshadeCenterDistance(d); }
void SiArrayConfig::SetLampshadeCenterAngleCmd(G4double angle) { SetLampshadeCenterAngle(angle); }
void SiArrayConfig::SetLampshadeHoleApothemCmd(G4double d) { SetLampshadeHoleApothem(d); }
void SiArrayConfig::SetLampshadeInnerDistanceCmd(G4double d) { SetLampshadeInnerDistance(d); }
void SiArrayConfig::SetLampshadeLengthCmd(G4double d) { SetLampshadeLength(d); }
void SiArrayConfig::SetLampshadeOuterAngleCmd(G4double angle) { SetLampshadeOuterAngle(angle); }
void SiArrayConfig::SetLampshadeModuleSpanCmd(G4double angle) { SetLampshadeModuleSpan(angle); }
void SiArrayConfig::SetLampshadeThicknessCmd(G4double d) { SetLampshadeThickness(d); }
void SiArrayConfig::SetCapFoldAngleCmd(G4double angle) { SetCapFoldAngle(angle); }
void SiArrayConfig::SetCapClearanceCmd(G4double d) { SetCapClearance(d); }
void SiArrayConfig::SetCapStaggerCmd(G4double d) { SetCapStagger(d); }
void SiArrayConfig::SetForwardAlphaTestBeamCmd(G4bool enabled) { SetForwardAlphaTestBeam(enabled); }
void SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd(G4double energy)
{
  SetForwardAlphaTestBeamEnergy(energy);
}
void SiArrayConfig::SetFixedAlphaTestBeamCmd(G4bool enabled) { SetFixedAlphaTestBeam(enabled); }
void SiArrayConfig::SetFixedAlphaTestBeamEnergyCmd(G4double energy)
{
  SetFixedAlphaTestBeamEnergy(energy);
}
void SiArrayConfig::SetFixedAlphaTestBeamThetaCmd(G4double theta) { SetFixedAlphaTestBeamTheta(theta); }
void SiArrayConfig::SetFixedAlphaTestBeamPhiCmd(G4double phi) { SetFixedAlphaTestBeamPhi(phi); }

void SiArrayConfig::PrintConfigCommand()
{
  const G4double s3_inner_angle = std::atan2(11. * mm, g_forward_distance);
  const G4double s3_outer_angle = std::atan2(35. * mm, g_forward_distance);
  const G4double surface_edge_distance = std::hypot(GetLampshadeInnerRadius(),
                                                    GetLampshadeInnerDistance());

  G4cout << "\n===== Si array runtime configuration (inclined lampshade + rear S3) =====" << G4endl
         << "  W1 drum modules                    = " << g_drum_modules << G4endl
         << "  W1 drum inscribed radius           = " << GetDrumInscribedRadius() / mm << " mm" << G4endl
         << "  W1 strips/module                   = " << g_barrel_strips_phi << " front + "
         << g_barrel_strips_z << " back" << G4endl
         << "  lampshade modules/end              = " << g_lampshade_modules << G4endl
         << "  lampshade strips/module            = " << g_lampshade_sectors << " sectors + "
         << g_lampshade_rings << " rings" << G4endl
         << "  lampshade module phi span          = " << g_lampshade_module_span / deg << " deg" << G4endl
         << "  lampshade centre distance/angle    = " << g_lampshade_center_distance / mm << " mm / "
         << g_lampshade_center_angle / deg << " deg" << G4endl
         << "  lampshade narrow edge (r,z,theta)  = " << GetLampshadeInnerRadius() / mm << ", "
         << GetLampshadeInnerDistance() / mm << " mm, "
         << GetLampshadeInnerAngle() / deg << " deg" << G4endl
         << "  lampshade wide edge (r,z,theta)    = " << GetLampshadeOuterRadius() / mm << ", "
         << GetLampshadeOuterDistance() / mm << " mm, "
         << GetLampshadeOuterAngle() / deg << " deg" << G4endl
         << "  lampshade edge target distance     = " << surface_edge_distance / mm << " mm" << G4endl
         << "  lampshade active length/thickness  = " << g_lampshade_length / mm << " / "
         << g_lampshade_thickness / mm << " mm" << G4endl
         << "  forward/backward S3 distance       = " << g_forward_distance / mm << " / "
         << g_backward_distance / mm << " mm" << G4endl
         << "  forward S3 angular coverage        = " << s3_inner_angle / deg << " - "
         << s3_outer_angle / deg << " deg" << G4endl
         << "  S3 strips/module                   = " << g_annular_sectors << " sectors + "
         << g_annular_rings << " rings" << G4endl
         << "  enableDrum                         = " << (g_enable_drum ? "true" : "false") << G4endl
         << "  enableForwardLampshade             = " << (g_enable_forward_cap ? "true" : "false") << G4endl
         << "  enableBackwardLampshade            = " << (g_enable_backward_cap ? "true" : "false") << G4endl
         << "  enableForwardS3                    = " << (g_enable_forward_annular ? "true" : "false") << G4endl
         << "  enableBackwardS3                   = " << (g_enable_backward_annular ? "true" : "false") << G4endl
         << "  fixedAlphaTestBeam                 = " << (g_fixed_alpha_test_beam ? "true" : "false")
         << " (" << g_fixed_alpha_test_beam_energy / MeV << " MeV, theta="
         << g_fixed_alpha_test_beam_theta / deg << " deg, phi="
         << g_fixed_alpha_test_beam_phi / deg << " deg)" << G4endl;
}

void SiArrayConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this,
                                                    "/si/",
                                                    "Si geometry and DSSD segmentation controls");

  messenger->DeclareMethod("barrelStripsZ", &SiArrayConfig::SetBarrelStripsZCmd,
                           "W1 back-face strips along local y");
  messenger->DeclareMethod("barrelStripsPhi", &SiArrayConfig::SetBarrelStripsPhiCmd,
                           "W1 front-face strips along local x");
  messenger->DeclareMethod("annularRings", &SiArrayConfig::SetAnnularRingsCmd,
                           "S3 radial rings");
  messenger->DeclareMethod("annularSectors", &SiArrayConfig::SetAnnularSectorsCmd,
                           "S3 angular sectors");
  messenger->DeclareMethod("lampshadeRings", &SiArrayConfig::SetLampshadeRingsCmd,
                           "Wedge-DSSD radial/ring strips");
  messenger->DeclareMethod("lampshadeSectors", &SiArrayConfig::SetLampshadeSectorsCmd,
                           "Wedge-DSSD transverse/sector strips");

  messenger->DeclareMethod("enableForwardAnnular", &SiArrayConfig::SetEnableForwardAnnularCmd,
                           "Enable S3 behind the forward central aperture");
  messenger->DeclareMethod("enableBackwardAnnular", &SiArrayConfig::SetEnableBackwardAnnularCmd,
                           "Enable S3 behind the backward central aperture");
  messenger->DeclareMethodWithUnit("forwardDistance", "mm", &SiArrayConfig::SetForwardDistanceCmd,
                                   "Forward S3 distance from target centre");
  messenger->DeclareMethodWithUnit("backwardDistance", "mm", &SiArrayConfig::SetBackwardDistanceCmd,
                                   "Backward S3 distance from target centre");

  messenger->DeclareMethod("drumModules", &SiArrayConfig::SetDrumModulesCmd,
                           "Number of W1 faces in the equatorial drum");
  messenger->DeclareMethod("enableDrum", &SiArrayConfig::SetEnableDrumCmd,
                           "Enable the equatorial W1 drum");

  // Keep historical command names for enable flags; their meaning is now
  // forward/backward wedge lampshade.
  messenger->DeclareMethod("enableForwardCap", &SiArrayConfig::SetEnableForwardCapCmd,
                           "Enable the forward twelve-wedge lampshade");
  messenger->DeclareMethod("enableBackwardCap", &SiArrayConfig::SetEnableBackwardCapCmd,
                           "Enable the backward twelve-wedge lampshade");
  messenger->DeclareMethod("lampshadeModules", &SiArrayConfig::SetLampshadeModulesCmd,
                           "Wedge DSSSDs per lampshade (default 12)");
  messenger->DeclareMethodWithUnit("lampshadeCenterDistance", "mm", &SiArrayConfig::SetLampshadeCenterDistanceCmd,
                                   "Target-to-centre distance of each wedge DSSSD");
  messenger->DeclareMethodWithUnit("lampshadeCenterAngle", "deg", &SiArrayConfig::SetLampshadeCenterAngleCmd,
                                   "Forward polar angle of the wedge centre; backward is mirrored");
  messenger->DeclareMethodWithUnit("lampshadeHoleApothem", "mm", &SiArrayConfig::SetLampshadeHoleApothemCmd,
                                   "Obsolete compatibility command; ignored");
  messenger->DeclareMethodWithUnit("lampshadeInnerDistance", "mm", &SiArrayConfig::SetLampshadeInnerDistanceCmd,
                                   "Obsolete compatibility command; ignored");
  messenger->DeclareMethodWithUnit("lampshadeLength", "mm", &SiArrayConfig::SetLampshadeLengthCmd,
                                   "Active wedge length in its local plane");
  messenger->DeclareMethodWithUnit("lampshadeOuterAngle", "deg", &SiArrayConfig::SetLampshadeOuterAngleCmd,
                                   "Obsolete compatibility command; ignored");
  messenger->DeclareMethodWithUnit("lampshadeModuleSpan", "deg", &SiArrayConfig::SetLampshadeModuleSpanCmd,
                                   "Azimuthal span of one wedge; must be below 360/modules");
  messenger->DeclareMethodWithUnit("lampshadeThickness", "mm", &SiArrayConfig::SetLampshadeThicknessCmd,
                                   "Wedge-DSSD silicon thickness");

  messenger->DeclareMethodWithUnit("capFoldAngle", "deg", &SiArrayConfig::SetCapFoldAngleCmd,
                                   "Obsolete compatibility command; ignored");
  messenger->DeclareMethodWithUnit("capClearance", "mm", &SiArrayConfig::SetCapClearanceCmd,
                                   "Obsolete compatibility command; ignored");
  messenger->DeclareMethodWithUnit("capStagger", "mm", &SiArrayConfig::SetCapStaggerCmd,
                                   "Obsolete compatibility command; ignored");

  messenger->DeclareMethod("forwardAlphaTestBeam", &SiArrayConfig::SetForwardAlphaTestBeamCmd,
                           "Detector test only: fire a forward alpha beam");
  messenger->DeclareMethodWithUnit("forwardAlphaTestBeamEnergy", "MeV", &SiArrayConfig::SetForwardAlphaTestBeamEnergyCmd,
                                   "Forward-alpha test energy");
  messenger->DeclareMethod("fixedAlphaTestBeam", &SiArrayConfig::SetFixedAlphaTestBeamCmd,
                           "Detector test only: one alpha from the target centre at fixed theta/phi");
  messenger->DeclareMethodWithUnit("fixedAlphaTestBeamEnergy", "MeV", &SiArrayConfig::SetFixedAlphaTestBeamEnergyCmd,
                                   "Fixed-alpha test energy");
  messenger->DeclareMethodWithUnit("fixedAlphaTestBeamTheta", "deg", &SiArrayConfig::SetFixedAlphaTestBeamThetaCmd,
                                   "Fixed-alpha lab polar angle");
  messenger->DeclareMethodWithUnit("fixedAlphaTestBeamPhi", "deg", &SiArrayConfig::SetFixedAlphaTestBeamPhiCmd,
                                   "Fixed-alpha lab azimuth");

  messenger->DeclareMethod("printConfig", &SiArrayConfig::PrintConfigCommand,
                           "Print current Si geometry and readout configuration");
}
