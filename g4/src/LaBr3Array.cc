#include "LaBr3Array.hh"
#include "LaBr3Detector.hh"

#include <cmath>
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

G4ThreeVector SidePortPosition(G4double radius, G4double z, G4int sector_id)
{
  switch (sector_id) {
  case 0:
    return G4ThreeVector(radius, 0., z);
  case 1:
    return G4ThreeVector(-radius, 0., z);
  case 2:
    return G4ThreeVector(0., radius, z);
  default:
    return G4ThreeVector(0., -radius, z);
  }
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Array::LaBr3Array(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  labr3_numbers = 0;

  for (auto it = LaBr3Detector::map_name_to_sectors.begin(); it != LaBr3Detector::map_name_to_sectors.end(); it++) {
    labr3_numbers += it->second;
  }

  //
  G4NistManager* nist_manager = G4NistManager::Instance();

  G4Element* la_element = nist_manager->FindOrBuildElement("La");
  G4Element* br_element = nist_manager->FindOrBuildElement("Br");

  labr3_mat = new G4Material("G4_LABR3", 5.08 * g / cm3, 2);
  labr3_mat->AddElement(la_element, 1);
  labr3_mat->AddElement(br_element, 3);

  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");

  PrintDetectorDimensionInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Array::~LaBr3Array() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Array::Construct()
{
  std::vector<LaBr3Detector*>::iterator it = v_labr3_detector.begin();
  for (; it != v_labr3_detector.end(); it++)
    delete *it;
  v_labr3_detector.clear();

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  G4int ii = 0;
  for (auto map_it = LaBr3Detector::map_name_to_ring_id.begin(); map_it != LaBr3Detector::map_name_to_ring_id.end(); map_it++) {
    for (G4int j = 0; j < LaBr3Detector::map_name_to_sectors[map_it->first]; j++) {
      map_i2ring[ii] = map_it->second;
      map_i2sector[ii] = j;
      map_i2name[ii] = map_it->first;
      ii++;
    }
  }

  for (auto map_it = map_i2ring.begin(); map_it != map_i2ring.end(); map_it++) {
    G4cout << "i " << map_it->first << " ring " << map_it->second << G4endl;
  }
  for (auto map_it = map_i2sector.begin(); map_it != map_i2sector.end(); map_it++) {
    G4cout << "i " << map_it->first << " sector " << map_it->second << G4endl;
  }
  for (auto map_it = map_i2name.begin(); map_it != map_i2name.end(); map_it++) {
    G4cout << "i " << map_it->first << " name " << map_it->second << G4endl;
  }

  for (G4int i = 0; i < labr3_numbers; i++) {
    v_labr3_detector.push_back(new LaBr3Detector(exp_hall_log));
    v_labr3_detector[i]->SetName(map_i2name[i]);
    v_labr3_detector[i]->SetRingId(map_i2ring[i]);
    v_labr3_detector[i]->SetSectorId(map_i2sector[i]);
  }

  for (it = v_labr3_detector.begin(); it != v_labr3_detector.end(); it++) {
    (*it)->ConstructLaBr3Detector(LaBr3Detector::map_labr3_par[(*it)->GetName()], labr3_mat);
    (*it)->PlaceLaBr3Detector(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
    (*it)->ConstructAlShell(LaBr3Detector::map_al_par[(*it)->GetName()], al_mat);
    (*it)->PlaceAlShell(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Array::MakeSensitive(LaBr3SD* labr3_sd)
{
  std::vector<LaBr3Detector*>::iterator it = v_labr3_detector.begin();
  for (; it != v_labr3_detector.end(); it++) {
    (*it)->GetLog()->SetSensitiveDetector(labr3_sd);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Transform3D LaBr3Array::CalculatePlacement(G4String name, G4int sector_id)
{
  const G4double radius = std::hypot(LaBr3Detector::map_placement_par[name][0], LaBr3Detector::map_placement_par[name][1]) * mm;
  const G4double z = LaBr3Detector::map_placement_par[name][2] * mm;
  const G4ThreeVector pos = SidePortPosition(radius, z, sector_id);
  const G4ThreeVector target_pos(0., 0., TargetZPos);
  return G4Transform3D(MakeRotationLocalZToDirection(target_pos - pos), pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Array::PrintDetectorDimensionInfo()
{
  for (auto it = LaBr3Detector::map_labr3_par.begin(); it != LaBr3Detector::map_labr3_par.end(); it++) {
    std::cout << it->first << std::endl;
    for (auto j = 0; j < it->second.size(); j++) {
      std::cout << it->second[j] << " ";
    }
    std::cout << std::endl;
  }
}
