#include "SiArray.hh"

#include "SiArrayConfig.hh"
#include "SiDetector.hh"

#include "G4NistManager.hh"
#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "TString.h"

#include <cmath>
#include <string>
#include <utility>

namespace {
G4RotationMatrix MakeRotationLocalZToDirection(const G4ThreeVector& direction)
{
  const G4ThreeVector local_z = direction.unit();
  const G4ThreeVector reference = std::abs(local_z.z()) < 0.95
                                      ? G4ThreeVector(0., 0., 1.)
                                      : G4ThreeVector(0., 1., 0.);
  const G4ThreeVector local_x = reference.cross(local_z).unit();
  const G4ThreeVector local_y = local_z.cross(local_x).unit();
  return G4RotationMatrix(local_x, local_y, local_z);
}

G4bool IsLampshade(const G4String& name)
{
  return name == "Si_ForwardCap" || name == "Si_BackwardCap";
}
} // namespace

std::vector<SiArray::DetectorGeometry> SiArray::detector_geometry;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiArray::SiArray(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  si_numbers = 0;
  for (const auto& item : SiDetector::map_name_to_sectors) {
    si_numbers += item.second;
  }

  auto* nist = G4NistManager::Instance();
  si_mat = nist->FindOrBuildMaterial("G4_Si");

  PrintDetectorDimensionInfo();
}

SiArray::~SiArray() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::Construct()
{
  for (auto* detector : v_si_detector) delete detector;
  v_si_detector.clear();
  detector_geometry.clear();

  // Refresh the runtime-derived wedge dimensions before constructing solids.
  const std::array<G4double, 4> wedge_parameters = {
      2. * SiArrayConfig::GetLampshadeInnerHalfWidth() / mm,
      2. * SiArrayConfig::GetLampshadeOuterHalfWidth() / mm,
      SiArrayConfig::GetLampshadeLength() / mm,
      SiArrayConfig::GetLampshadeThickness() / mm};
  SiDetector::map_si_wedge_par["Si_ForwardCap"] = wedge_parameters;
  SiDetector::map_si_wedge_par["Si_BackwardCap"] = wedge_parameters;
  SiDetector::map_name_to_sectors["Si_Drum"] = SiArrayConfig::GetDrumModules();
  SiDetector::map_name_to_sectors["Si_ForwardCap"] = SiArrayConfig::GetLampshadeModules();
  SiDetector::map_name_to_sectors["Si_BackwardCap"] = SiArrayConfig::GetLampshadeModules();

  std::map<G4int, G4int> index_to_subarray;
  std::map<G4int, G4int> index_to_module;
  std::map<G4int, G4String> index_to_name;

  G4int index = 0;
  for (const auto& item : SiDetector::map_name_to_ring_id) {
    const G4String& name = item.first;
    if (name == "Si_ForwardAnnular" && !SiArrayConfig::GetEnableForwardAnnular()) continue;
    if (name == "Si_BackwardAnnular" && !SiArrayConfig::GetEnableBackwardAnnular()) continue;
    if (name == "Si_Drum" && !SiArrayConfig::GetEnableDrum()) continue;
    if (name == "Si_ForwardCap" && !SiArrayConfig::GetEnableForwardCap()) continue;
    if (name == "Si_BackwardCap" && !SiArrayConfig::GetEnableBackwardCap()) continue;

    G4int modules = SiDetector::map_name_to_sectors[name];
    if (name == "Si_Drum") modules = SiArrayConfig::GetDrumModules();
    if (IsLampshade(name)) modules = SiArrayConfig::GetLampshadeModules();

    for (G4int module = 0; module < modules; ++module) {
      index_to_subarray[index] = item.second;
      index_to_module[index] = module;
      index_to_name[index] = name;
      ++index;
    }
  }
  si_numbers = index;

  for (G4int i = 0; i < si_numbers; ++i) {
    auto* detector = new SiDetector(exp_hall_log);
    detector->SetName(index_to_name[i]);
    detector->SetRingId(index_to_subarray[i]);
    detector->SetSectorId(index_to_module[i]);
    v_si_detector.push_back(detector);
  }

  for (auto* detector : v_si_detector) {
    const G4String name = detector->GetName();

    const auto box_it = SiDetector::map_si_box_par.find(name);
    const auto wedge_it = SiDetector::map_si_wedge_par.find(name);
    const auto annular_it = SiDetector::map_si_par.find(name);

    if (box_it != SiDetector::map_si_box_par.end()) {
      detector->ConstructSiBoxDetector(box_it->second, si_mat);
    } else if (wedge_it != SiDetector::map_si_wedge_par.end()) {
      detector->ConstructSiWedgeDetector(wedge_it->second, si_mat);
    } else if (annular_it != SiDetector::map_si_par.end()) {
      detector->ConstructSiDetector(annular_it->second, si_mat);
    } else {
      G4ExceptionDescription description;
      description << "No silicon solid parameters are registered for " << name;
      G4Exception("SiArray::Construct", "SiArray001", FatalException, description);
    }

    const G4Transform3D placement = CalculatePlacement(name, detector->GetSectorId());
    detector->PlaceSiDetector(placement);
    RegisterDetectorGeometry(*detector, placement);

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::MakeSensitive(SiSD* si_sd)
{
  for (auto* detector : v_si_detector) {
    detector->GetLog()->SetSensitiveDetector(si_sd);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Transform3D SiArray::CalculatePlacement(G4String name, G4int module_id)
{
  const G4ThreeVector target_position(0., 0., TargetZPos);

  if (name == "Si_Drum") {
    const G4int modules = SiArrayConfig::GetDrumModules();
    const G4double radius = SiArrayConfig::GetDrumInscribedRadius();
    const G4double phi = twopi * static_cast<G4double>(module_id)
                         / static_cast<G4double>(modules);
    const G4ThreeVector position(radius * std::cos(phi),
                                 radius * std::sin(phi),
                                 TargetZPos);

    // Local z faces the target; local x is approximately azimuthal and local
    // y approximately follows the beam axis.
    return G4Transform3D(MakeRotationLocalZToDirection(target_position - position),
                         position);
  }

  if (IsLampshade(name)) {
    const G4bool forward = name == "Si_ForwardCap";
    const G4int modules = SiArrayConfig::GetLampshadeModules();
    const G4double phi = twopi * static_cast<G4double>(module_id)
                         / static_cast<G4double>(modules);
    const G4double theta_forward = SiArrayConfig::GetLampshadeCenterAngle();
    const G4double theta = forward ? theta_forward : pi - theta_forward;
    const G4double centre_distance = SiArrayConfig::GetLampshadeCenterDistance();

    const G4ThreeVector radial_out(std::cos(phi), std::sin(phi), 0.);
    const G4ThreeVector beam_axis(0., 0., 1.);
    const G4ThreeVector spherical_radial = std::sin(theta) * radial_out
                                           + std::cos(theta) * beam_axis;
    const G4ThreeVector theta_direction = std::cos(theta) * radial_out
                                          - std::sin(theta) * beam_axis;
    const G4ThreeVector position = target_position
                                   + centre_distance * spherical_radial;

    // Narrow edge is local y=-L/2 and must face the beam axis on both ends.
    // Forward: local +y follows increasing theta. Backward: local +y follows
    // decreasing theta, so the narrow edge remains near theta=180 degrees.
    const G4ThreeVector local_y = forward ? theta_direction : -theta_direction;
    G4ThreeVector local_x(-std::sin(phi), std::cos(phi), 0.);
    G4ThreeVector local_z = local_x.cross(local_y).unit();
    if (local_z.dot(target_position - position) < 0.) {
      local_x = -local_x;
      local_z = -local_z;
    }

    return G4Transform3D(G4RotationMatrix(local_x, local_y, local_z), position);
  }

  const G4double x = SiDetector::map_placement_par[name][0] * mm;
  const G4double y = SiDetector::map_placement_par[name][1] * mm;
  G4double z = SiDetector::map_placement_par[name][2] * mm;

  if (name == "Si_BackwardAnnular") {
    z = TargetZPos - SiArrayConfig::GetBackwardDistance();
  } else if (name == "Si_ForwardAnnular") {
    z = TargetZPos + SiArrayConfig::GetForwardDistance();
  }

  const G4ThreeVector position(x, y, z);
  return G4Transform3D(MakeRotationLocalZToDirection(target_position - position),
                       position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::RegisterDetectorGeometry(SiDetector& detector,
                                       const G4Transform3D& transform)
{
  DetectorGeometry geometry;
  geometry.detector_id = detector.GetPhy()->GetCopyNo();
  geometry.subarray_id = detector.GetRingId();
  geometry.module_id = detector.GetSectorId();
  geometry.local_to_world = transform;

  const G4String name = detector.GetName();
  const auto box_it = SiDetector::map_si_box_par.find(name);
  const auto wedge_it = SiDetector::map_si_wedge_par.find(name);
  const auto annular_it = SiDetector::map_si_par.find(name);

  if (box_it != SiDetector::map_si_box_par.end()) {
    geometry.detector_model = 1;
    geometry.half_x = box_it->second[0] / 2. * mm;
    geometry.half_y = box_it->second[1] / 2. * mm;
  } else if (annular_it != SiDetector::map_si_par.end()) {
    geometry.detector_model = 2;
    geometry.r_outer = annular_it->second[0] / 2. * mm;
    const auto inner_it = SiDetector::map_si_inner_radius.find(name);
    if (inner_it != SiDetector::map_si_inner_radius.end()) {
      geometry.r_inner = inner_it->second * mm;
    }
  } else if (wedge_it != SiDetector::map_si_wedge_par.end()) {
    geometry.detector_model = 3;
    geometry.wedge_inner_half_width = wedge_it->second[0] / 2. * mm;
    geometry.wedge_outer_half_width = wedge_it->second[1] / 2. * mm;
    geometry.wedge_half_length = wedge_it->second[2] / 2. * mm;
  }

  detector_geometry.push_back(geometry);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<SiPixelMapData> SiArray::BuildPixelMap()
{
  std::vector<SiPixelMapData> pixel_map;

  for (const auto& geometry : detector_geometry) {
    G4int front_count = 0;
    G4int back_count = 0;

    if (geometry.detector_model == 1) {
      front_count = SiArrayConfig::GetBarrelStripsPhi();
      back_count = SiArrayConfig::GetBarrelStripsZ();
    } else if (geometry.detector_model == 2) {
      front_count = SiArrayConfig::GetAnnularSectors();
      back_count = SiArrayConfig::GetAnnularRings();
    } else if (geometry.detector_model == 3) {
      front_count = SiArrayConfig::GetLampshadeSectors();
      back_count = SiArrayConfig::GetLampshadeRings();
    }

    if (front_count <= 0 || back_count <= 0) continue;
    pixel_map.reserve(pixel_map.size()
                      + static_cast<std::size_t>(front_count * back_count));

    for (G4int front_strip = 0; front_strip < front_count; ++front_strip) {
      for (G4int back_strip = 0; back_strip < back_count; ++back_strip) {
        G4ThreeVector local_center;

        if (geometry.detector_model == 1) {
          const G4double x = -geometry.half_x
              + (static_cast<G4double>(front_strip) + 0.5)
                    * (2. * geometry.half_x / front_count);
          const G4double y = -geometry.half_y
              + (static_cast<G4double>(back_strip) + 0.5)
                    * (2. * geometry.half_y / back_count);
          local_center = G4ThreeVector(x, y, 0.);
        } else if (geometry.detector_model == 2) {
          const G4double phi = (static_cast<G4double>(front_strip) + 0.5)
                               * twopi / front_count;
          const G4double radius = geometry.r_inner
              + (static_cast<G4double>(back_strip) + 0.5)
                    * ((geometry.r_outer - geometry.r_inner) / back_count);
          local_center = G4ThreeVector(radius * std::cos(phi),
                                       radius * std::sin(phi),
                                       0.);
        } else {
          const G4double y = -geometry.wedge_half_length
              + (static_cast<G4double>(back_strip) + 0.5)
                    * (2. * geometry.wedge_half_length / back_count);
          const G4double half_width = SiArrayConfig::LampshadeHalfWidthAtY(
              y,
              geometry.wedge_half_length,
              geometry.wedge_inner_half_width,
              geometry.wedge_outer_half_width);
          const G4double x = -half_width
              + (static_cast<G4double>(front_strip) + 0.5)
                    * (2. * half_width / front_count);
          local_center = G4ThreeVector(x, y, 0.);
        }

        const G4ThreeVector world_center =
            geometry.local_to_world.getRotation() * local_center
            + geometry.local_to_world.getTranslation();

        const G4double dx = world_center.x();
        const G4double dy = world_center.y();
        const G4double dz = world_center.z() - TargetZPos;
        const G4double theta = std::atan2(std::hypot(dx, dy), dz) / deg;
        G4double phi = std::atan2(dy, dx) / deg;
        if (phi < 0.) phi += 360.;

        SiPixelMapData data;
        data.detector_model = geometry.detector_model;
        data.subarray_id = geometry.subarray_id;
        data.module_id = geometry.module_id;
        data.front_strip_id = front_strip;
        data.back_strip_id = back_strip;
        data.x_center_mm = world_center.x() / mm;
        data.y_center_mm = world_center.y() / mm;
        data.z_center_mm = world_center.z() / mm;
        data.theta_lab_center_deg = theta;
        data.phi_lab_center_deg = phi;
        pixel_map.push_back(data);
      }
    }
  }

  return pixel_map;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::PrintDetectorDimensionInfo()
{
  G4cout << "[SiArray] annular detector parameters" << G4endl;
  for (const auto& item : SiDetector::map_si_par) {
    G4cout << "  " << item.first << ": outer diameter=" << item.second[0]
           << " mm, thickness=" << item.second[1] << " mm" << G4endl;
  }
  G4cout << "[SiArray] wedge lampshade: narrow width="
         << 2. * SiArrayConfig::GetLampshadeInnerHalfWidth() / mm
         << " mm, wide width="
         << 2. * SiArrayConfig::GetLampshadeOuterHalfWidth() / mm
         << " mm, length=" << SiArrayConfig::GetLampshadeLength() / mm
         << " mm" << G4endl;
}
