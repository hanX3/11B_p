#include "LaBr3Array.hh"
#include "LaBr3Detector.hh"

#include <cmath>
#include "TString.h"

//
LaBr3Array::LaBr3Array(G4LogicalVolume *log)
: exp_hall_log(log)
{
  labr3_numbers = 0;

  for(auto it=LaBr3Detector::map_name_to_sectors.begin();it!=LaBr3Detector::map_name_to_sectors.end();it++){
    labr3_numbers += it->second;
  }
  
  //
  G4NistManager *nist_manager = G4NistManager::Instance();

  G4Element *la_element = nist_manager->FindOrBuildElement("La");
  G4Element *br_element = nist_manager->FindOrBuildElement("Br");

  labr3_mat = new G4Material("G4_LABR3", 5.08 *g/cm3, 2);
  labr3_mat->AddElement(la_element, 1);
  labr3_mat->AddElement(br_element, 3);

  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");
  
  PrintDetectorDimensionInfo();
}

//
LaBr3Array::~LaBr3Array()
{

}

//
void LaBr3Array::Construct()
{
  std::vector<LaBr3Detector*>::iterator it = v_labr3_detector.begin();
  // clear all elements from the array
  for(;it!=v_labr3_detector.end();it++) delete *it;
  v_labr3_detector.clear();

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  int ii = 0;
  for(auto it=LaBr3Detector::map_name_to_ring_id.begin();it!=LaBr3Detector::map_name_to_ring_id.end();it++){
    for(int j=0;j<LaBr3Detector::map_name_to_sectors[it->first];j++){
      map_i2ring[ii] = it->second;
      map_i2sector[ii] = j;
      map_i2name[ii] = it->first;
      ii++;
    }
  }

  for(auto it=map_i2ring.begin();it!=map_i2ring.end();it++){
    G4cout << "i " << it->first << " ring " << it->second << G4endl;
  }
  for(auto it=map_i2sector.begin();it!=map_i2sector.end();it++){
    G4cout << "i " << it->first << " sector " << it->second << G4endl;
  }
  for(auto it=map_i2name.begin();it!=map_i2name.end();it++){
    G4cout << "i " << it->first << " name " << it->second << G4endl;
  }

  for(int i=0;i<labr3_numbers;i++){
    v_labr3_detector.push_back(new LaBr3Detector(exp_hall_log));
    v_labr3_detector[i]->SetName(map_i2name[i]);
    v_labr3_detector[i]->SetRingId(map_i2ring[i]);
    v_labr3_detector[i]->SetSectorId(map_i2sector[i]);
  }

  for(it=v_labr3_detector.begin();it!=v_labr3_detector.end();it++){
    (*it)->ConstructLaBr3Detector(LaBr3Detector::map_labr3_par[(*it)->GetName()], labr3_mat);
    (*it)->PlaceLaBr3Detector(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));

    (*it)->ConstructAlShell(LaBr3Detector::map_al_par[(*it)->GetName()], al_mat);
    (*it)->PlaceAlShell(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
  }
}

//
void LaBr3Array::MakeSensitive(LaBr3SD *labr3_sd)
{
  std::vector<LaBr3Detector*>::iterator it = v_labr3_detector.begin();
  for(;it!=v_labr3_detector.end();it++){
     (*it)->GetLog()->SetSensitiveDetector(labr3_sd);
  }
}

//
G4Transform3D LaBr3Array::CalculatePlacement(G4String name, G4int sector_id)
{
  G4RotationMatrix *rot_matrix = new G4RotationMatrix();
  G4double rot_x_angle = 0. *deg;
  G4double rot_y_angle = 90. *deg;
  G4double rot_z_angle = 0. *deg;
  G4cout << "x angle " << rot_x_angle << " y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  rot_matrix->rotateX(rot_x_angle);
  rot_matrix->rotateY(rot_y_angle);
  rot_matrix->rotateZ(rot_z_angle);

  rot_matrix->print(G4cout);

  G4double x = LaBr3Detector::map_placement_par[name][0] *mm;
  G4double y = LaBr3Detector::map_placement_par[name][1] *mm;
  G4double z = LaBr3Detector::map_placement_par[name][2] *mm;
  G4ThreeVector pos = G4ThreeVector(x, y, z);
  G4cout << "pos x " << pos.x() << "pos y " << pos.y() << "pos z " << pos.z() << G4endl;

  G4Transform3D transform(*rot_matrix, pos);

  return transform;
}

//
void LaBr3Array::PrintDetectorDimensionInfo()
{
  for(auto it=LaBr3Detector::map_labr3_par.begin();it!=LaBr3Detector::map_labr3_par.end();it++){
    std::cout << it->first << std::endl;
    for(auto j=0;j<it->second.size();j++){
      std::cout << it->second[j] << " ";
    }
    std::cout <<std::endl;
  }
}