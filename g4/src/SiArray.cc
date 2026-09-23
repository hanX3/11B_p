#include "SiArray.hh"
#include "SiDetector.hh"
#include "SiArrayConfig.hh"

#include <cmath>
#include "G4PhysicalConstants.hh"
#include "TString.h"

namespace {
G4RotationMatrix MakeBarrelModuleRotation(G4double phi)
{
  const G4ThreeVector local_y(0., 0., 1.);
  const G4ThreeVector local_z(-std::cos(phi), -std::sin(phi), 0.);
  const G4ThreeVector local_x = local_y.cross(local_z).unit();
  return G4RotationMatrix(local_x, local_y, local_z);
}

G4RotationMatrix MakeRotationLocalZToDirection(const G4ThreeVector& direction)
{
  const G4ThreeVector local_z = direction.unit();
  const G4ThreeVector reference = std::abs(local_z.z()) < 0.95 ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(0., 1., 0.);
  const G4ThreeVector local_x = reference.cross(local_z).unit();
  const G4ThreeVector local_y = local_z.cross(local_x).unit();
  return G4RotationMatrix(local_x, local_y, local_z);
}
} // namespace

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

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  G4int ii = 0;
  for (auto map_it = SiDetector::map_name_to_ring_id.begin(); map_it != SiDetector::map_name_to_ring_id.end(); map_it++) {
    // The forward annular DSSD is optional; skip building it when disabled so
    // the geometry and sensitive-detector list stay consistent.
    if (map_it->first == "Si_ForwardAnnular" && !SiArrayConfig::GetEnableForwardAnnular()) continue;
    for (G4int j = 0; j < SiDetector::map_name_to_sectors[map_it->first]; j++) {
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
    if ((*it)->GetName() == "Si_Barrel") {
      (*it)->ConstructSiBoxDetector(SiDetector::map_si_box_par[(*it)->GetName()], si_mat);
    } else {
      (*it)->ConstructSiDetector(SiDetector::map_si_par[(*it)->GetName()], si_mat);
    }
    (*it)->PlaceSiDetector(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
    if (SiDetector::map_al_par.find((*it)->GetName()) != SiDetector::map_al_par.end()) {
      (*it)->ConstructAlShell(SiDetector::map_al_par[(*it)->GetName()], al_mat);
      (*it)->PlaceAlShell(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
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

  if (name == "Si_Barrel") {
    const G4double radius = SiDetector::map_placement_par[name][0] * mm;
    const G4int n_modules = SiDetector::map_name_to_sectors[name];
    const G4double phi = twopi * static_cast<G4double>(sector_id) / static_cast<G4double>(n_modules);
    const G4ThreeVector pos(radius * std::cos(phi), radius * std::sin(phi), TargetZPos);
    return G4Transform3D(MakeBarrelModuleRotation(phi), pos);
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
