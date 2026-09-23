#include "SiArray.hh"
#include "SiDetector.hh"
#include "SiArrayConfig.hh"

#include <cmath>
#include <utility>
#include "G4PhysicalConstants.hh"
#include "TString.h"

namespace {
G4RotationMatrix MakeRotationLocalZToDirection(const G4ThreeVector& direction)
{
  const G4ThreeVector local_z = direction.unit();
  const G4ThreeVector reference = std::abs(local_z.z()) < 0.95 ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(0., 1., 0.);
  const G4ThreeVector local_x = reference.cross(local_z).unit();
  const G4ThreeVector local_y = local_z.cross(local_x).unit();
  return G4RotationMatrix(local_x, local_y, local_z);
}
} // namespace

std::vector<SiArray::DetectorGeometry> SiArray::detector_geometry;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiArray::SiArray(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  si_numbers = 0;

  for (auto it = SiDetector::map_name_to_sectors.begin(); it != SiDetector::map_name_to_sectors.end(); it++) {
    si_numbers += it->second;
  }

  //
  G4NistManager* nist_manager = G4NistManager::Instance();

  si_mat = nist_manager->FindOrBuildMaterial("G4_Si");
  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");

  PrintDetectorDimensionInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiArray::~SiArray() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::Construct()
{
  std::vector<SiDetector*>::iterator it = v_si_detector.begin();
  for (; it != v_si_detector.end(); it++)
    delete *it;
  v_si_detector.clear();
  detector_geometry.clear();

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  G4int ii = 0;
  for (auto map_it = SiDetector::map_name_to_ring_id.begin(); map_it != SiDetector::map_name_to_ring_id.end(); map_it++) {
    // Optional sub-arrays: skip disabled ones so the geometry and
    // sensitive-detector list stay consistent.
    if (map_it->first == "Si_ForwardAnnular" && !SiArrayConfig::GetEnableForwardAnnular()) continue;
    if (map_it->first == "Si_Drum" && !SiArrayConfig::GetEnableDrum()) continue;
    if (map_it->first == "Si_ForwardCap" && !SiArrayConfig::GetEnableForwardCap()) continue;
    if (map_it->first == "Si_BackwardCap" && !SiArrayConfig::GetEnableBackwardCap()) continue;
    // Drum and cap module counts are runtime-configurable (caps always match
    // the drum face count so the hinges line up); annulars are single modules.
    const G4bool is_w1_ring = map_it->first == "Si_Drum" || map_it->first.find("Cap") != std::string::npos;
    const G4int n_modules = is_w1_ring ? SiArrayConfig::GetDrumModules() : SiDetector::map_name_to_sectors[map_it->first];
    for (G4int j = 0; j < n_modules; j++) {
      map_i2ring[ii] = map_it->second;
      map_i2sector[ii] = j;
      map_i2name[ii] = map_it->first;
      ii++;
    }
  }
  si_numbers = ii;

  for (auto map_it = map_i2ring.begin(); map_it != map_i2ring.end(); map_it++) {
    G4cout << "i " << map_it->first << " ring " << map_it->second << G4endl;
  }
  for (auto map_it = map_i2sector.begin(); map_it != map_i2sector.end(); map_it++) {
    G4cout << "i " << map_it->first << " sector " << map_it->second << G4endl;
  }
  for (auto map_it = map_i2name.begin(); map_it != map_i2name.end(); map_it++) {
    G4cout << "i " << map_it->first << " name " << map_it->second << G4endl;
  }

  for (G4int i = 0; i < si_numbers; i++) {
    v_si_detector.push_back(new SiDetector(exp_hall_log));
    v_si_detector[i]->SetName(map_i2name[i]);
    v_si_detector[i]->SetRingId(map_i2ring[i]);
    v_si_detector[i]->SetSectorId(map_i2sector[i]);
  }

  for (it = v_si_detector.begin(); it != v_si_detector.end(); it++) {
    // Flat (box) W1 modules -- drum and caps -- vs annular (tubs) S3 DSSDs.
    if (SiDetector::map_si_box_par.find((*it)->GetName()) != SiDetector::map_si_box_par.end()) {
      (*it)->ConstructSiBoxDetector(SiDetector::map_si_box_par[(*it)->GetName()], si_mat);
    } else {
      (*it)->ConstructSiDetector(SiDetector::map_si_par[(*it)->GetName()], si_mat);
    }

    const G4Transform3D placement = CalculatePlacement((*it)->GetName(), (*it)->GetSectorId());
    (*it)->PlaceSiDetector(placement);
    RegisterDetectorGeometry(*(*it), placement);

    if (SiDetector::map_al_par.find((*it)->GetName()) != SiDetector::map_al_par.end()) {
      (*it)->ConstructAlShell(SiDetector::map_al_par[(*it)->GetName()], al_mat);
      (*it)->PlaceAlShell(placement);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::MakeSensitive(SiSD* si_sd)
{
  std::vector<SiDetector*>::iterator it = v_si_detector.begin();
  for (; it != v_si_detector.end(); it++) {
    (*it)->GetLog()->SetSensitiveDetector(si_sd);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Transform3D SiArray::CalculatePlacement(G4String name, G4int sector_id)
{
  const G4ThreeVector target_pos(0., 0., TargetZPos);
  const G4double half_width = 25. * mm; // W1 active half-width

  if (name == "Si_Drum") {
    // Regular n-gon prism: face centers on the inscribed circle, face normals
    // radial.  Faces close exactly because the inscribed radius is derived
    // from the module count (SiArrayConfig::GetDrumInscribedRadius).
    const G4int n = SiArrayConfig::GetDrumModules();
    const G4double r_in = SiArrayConfig::GetDrumInscribedRadius();
    const G4double phi = twopi * static_cast<G4double>(sector_id) / static_cast<G4double>(n);
    const G4ThreeVector pos(r_in * std::cos(phi), r_in * std::sin(phi), TargetZPos);
    // Local frame: z -> inward radial, x -> azimuthal (phi strips),
    // y -> beam axis (polar strips).
    return G4Transform3D(MakeRotationLocalZToDirection(G4ThreeVector(0., 0., TargetZPos) - pos), pos);
  }

  if (name == "Si_ForwardCap" || name == "Si_BackwardCap") {
    // Fish-scale end cap.  Each square W1 plate is hinged on a drum end edge
    // and folded inward by capFoldAngle.  Folding full-width square plates on
    // EVERY face would make neighbours intersect, so the plates alternate
    // between two sub-rings:
    //   even sector_id -> inner sub-ring, hinge exactly on the drum edge
    //   odd  sector_id -> outer sub-ring, hinge capStagger further out and
    //                     (automatically, via sector phi) half a face period
    //                     rotated, overlapping the inner plates like scales.
    // A small clearance keeps the inner hinge from sharing a surface with the
    // drum solid (avoids touching-surface warnings from the overlap checker).
    const G4int n = SiArrayConfig::GetDrumModules();
    const G4double fold = SiArrayConfig::GetCapFoldAngle();
    const G4double stagger = SiArrayConfig::GetCapStagger();
    const G4double clearance = 0.3 * mm;
    const G4bool forward = (name == "Si_ForwardCap");
    const G4double axial_sign = forward ? +1. : -1.;

    const G4double phi = twopi * static_cast<G4double>(sector_id) / static_cast<G4double>(n);
    const G4double hinge_r = SiArrayConfig::GetDrumInscribedRadius() + ((sector_id % 2 == 1) ? stagger : 0.);
    const G4double hinge_z = TargetZPos + axial_sign * half_width;

    const G4ThreeVector radial_out(std::cos(phi), std::sin(phi), 0.);
    const G4ThreeVector tangent(-std::sin(phi), std::cos(phi), 0.);
    // "Up the plate" direction: inward in radius, away from the target axially.
    const G4ThreeVector up = (-std::sin(fold)) * radial_out + (axial_sign * std::cos(fold)) * G4ThreeVector(0., 0., 1.);

    const G4ThreeVector hinge_mid = hinge_r * radial_out + G4ThreeVector(0., 0., hinge_z);
    const G4ThreeVector pos = hinge_mid + (half_width + clearance) * up;

    // Local frame: x -> tangent (azimuthal strips), y -> up (polar strips),
    // z -> plate normal.  Same x/y convention as the drum, so
    // the W1 front/back strip-coordinate mapping applies unchanged.
    const G4ThreeVector local_x = tangent;
    const G4ThreeVector local_y = up.unit();
    const G4ThreeVector local_z = local_x.cross(local_y).unit();
    return G4Transform3D(G4RotationMatrix(local_x, local_y, local_z), pos);
  }

  const G4double x = SiDetector::map_placement_par[name][0] * mm;
  const G4double y = SiDetector::map_placement_par[name][1] * mm;
  // The axial distance of the two annular DSSDs is driven by SiArrayConfig so
  // it can be scanned from a macro; other detectors fall back to the map value.
  G4double z = SiDetector::map_placement_par[name][2] * mm;
  if (name == "Si_BackwardAnnular") {
    z = TargetZPos - SiArrayConfig::GetBackwardDistance();
  } else if (name == "Si_ForwardAnnular") {
    z = TargetZPos + SiArrayConfig::GetForwardDistance();
  }
  const G4ThreeVector pos(x, y, z);
  return G4Transform3D(MakeRotationLocalZToDirection(target_pos - pos), pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiArray::RegisterDetectorGeometry(SiDetector& detector, const G4Transform3D& transform)
{
  DetectorGeometry geometry;
  geometry.detector_id = detector.GetPhy()->GetCopyNo();
  geometry.subarray_id = detector.GetRingId();
  geometry.module_id = detector.GetSectorId();
  geometry.local_to_world = transform;

  const G4String name = detector.GetName();
  const auto box_it = SiDetector::map_si_box_par.find(name);
  if (box_it != SiDetector::map_si_box_par.end()) {
    geometry.detector_model = 1;
    geometry.half_x = box_it->second[0] / 2. * mm;
    geometry.half_y = box_it->second[1] / 2. * mm;
  } else {
    geometry.detector_model = 2;
    const auto tub_it = SiDetector::map_si_par.find(name);
    if (tub_it != SiDetector::map_si_par.end()) geometry.r_outer = tub_it->second[0] / 2. * mm;
    const auto inner_it = SiDetector::map_si_inner_radius.find(name);
    if (inner_it != SiDetector::map_si_inner_radius.end()) geometry.r_inner = inner_it->second * mm;
  }

  detector_geometry.push_back(geometry);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<SiPixelMapData> SiArray::BuildPixelMap()
{
  std::vector<SiPixelMapData> pixel_map;

  for (const auto& geometry : detector_geometry) {
    const G4int n_front = geometry.detector_model == 1 ? SiArrayConfig::GetBarrelStripsPhi()
                                                       : SiArrayConfig::GetAnnularSectors();
    const G4int n_back = geometry.detector_model == 1 ? SiArrayConfig::GetBarrelStripsZ()
                                                      : SiArrayConfig::GetAnnularRings();
    if (n_front <= 0 || n_back <= 0) continue;

    pixel_map.reserve(pixel_map.size() + static_cast<std::size_t>(n_front * n_back));

    for (G4int front_strip = 0; front_strip < n_front; ++front_strip) {
      for (G4int back_strip = 0; back_strip < n_back; ++back_strip) {
        G4ThreeVector local_center;
        if (geometry.detector_model == 1) {
          const G4double x = -geometry.half_x + (static_cast<G4double>(front_strip) + 0.5) *
                                                        (2. * geometry.half_x / n_front);
          const G4double y = -geometry.half_y + (static_cast<G4double>(back_strip) + 0.5) *
                                                        (2. * geometry.half_y / n_back);
          local_center = G4ThreeVector(x, y, 0.);
        } else {
          const G4double phi = (static_cast<G4double>(front_strip) + 0.5) * twopi / n_front;
          const G4double radius = geometry.r_inner + (static_cast<G4double>(back_strip) + 0.5) *
                                                               ((geometry.r_outer - geometry.r_inner) / n_back);
          local_center = G4ThreeVector(radius * std::cos(phi), radius * std::sin(phi), 0.);
        }

        const G4ThreeVector world_center = geometry.local_to_world.getRotation() * local_center +
                                           geometry.local_to_world.getTranslation();
        const G4double dx = world_center.x();
        const G4double dy = world_center.y();
        const G4double dz = world_center.z() - TargetZPos;
        const G4double theta = std::atan2(std::hypot(dx, dy), dz) / deg;
        G4double phi = std::atan2(dy, dx) / deg;
        if (phi < 0.) phi += 360.;

        SiPixelMapData data;
        data.detector_id = geometry.detector_id;
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
  for (auto it = SiDetector::map_si_par.begin(); it != SiDetector::map_si_par.end(); it++) {
    std::cout << it->first << std::endl;
    for (auto j = 0; j < it->second.size(); j++) {
      std::cout << it->second[j] << " ";
    }
    std::cout << std::endl;
  }
}
