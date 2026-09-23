#include "HPGeArray.hh"
#include "HPGeDetector.hh"

#include <cmath>
#include "TString.h"

//
HPGeArray::HPGeArray(G4LogicalVolume *log)
: exp_hall_log(log)
{
  hpge_numbers = 0;

  for(auto it=HPGeDetector::map_name_to_sectors.begin();it!=HPGeDetector::map_name_to_sectors.end();it++){
    hpge_numbers += it->second;
  }
  
  //
  G4NistManager *nist_manager = G4NistManager::Instance();

  hpge_mat = nist_manager->FindOrBuildMaterial("G4_Ge");
  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");
  
  PrintDetectorDimensionInfo();
}

//
HPGeArray::~HPGeArray()
{

}

//
void HPGeArray::Construct()
{
  std::vector<HPGeDetector*>::iterator it = v_hpge_detector.begin();
  // clear all elements from the array
  for(;it!=v_hpge_detector.end();it++) delete *it;
  v_hpge_detector.clear();

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  int ii = 0;
  for(auto it=HPGeDetector::map_name_to_ring_id.begin();it!=HPGeDetector::map_name_to_ring_id.end();it++){
    for(int j=0;j<HPGeDetector::map_name_to_sectors[it->first];j++){
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

  for(int i=0;i<hpge_numbers;i++){
    v_hpge_detector.push_back(new HPGeDetector(exp_hall_log));
    v_hpge_detector[i]->SetName(map_i2name[i]);
    v_hpge_detector[i]->SetRingId(map_i2ring[i]);
    v_hpge_detector[i]->SetSectorId(map_i2sector[i]);
  }

  for(auto it=v_hpge_detector.begin();it!=v_hpge_detector.end();it++){
    (*it)->ConstructHPGeDetector(HPGeDetector::map_hpge_par[(*it)->GetName()], hpge_mat);
    (*it)->PlaceHPGeDetector(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));

    (*it)->ConstructAlShell(HPGeDetector::map_al_par[(*it)->GetName()], al_mat);
    (*it)->PlaceAlShell(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
  }
}

//
void HPGeArray::MakeSensitive(HPGeSD *hpge_sd)
{
  std::vector<HPGeDetector*>::iterator it = v_hpge_detector.begin();
  for(;it!=v_hpge_detector.end();it++){
     (*it)->GetLog()->SetSensitiveDetector(hpge_sd);
  }
}

//
G4Transform3D HPGeArray::CalculatePlacement(G4String name, G4int sector_id)
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

  G4double x = HPGeDetector::map_placement_par[name][0] *mm;
  G4double y = HPGeDetector::map_placement_par[name][1] *mm;
  G4double z = HPGeDetector::map_placement_par[name][2] *mm;
  G4ThreeVector pos = G4ThreeVector(x, y, z);
  G4cout << "pos x " << pos.x() << "pos y " << pos.y() << "pos z " << pos.z() << G4endl;

  G4Transform3D transform(*rot_matrix, pos);

  return transform;
}

//
void HPGeArray::PrintDetectorDimensionInfo()
{
  for(auto it=HPGeDetector::map_hpge_par.begin();it!=HPGeDetector::map_hpge_par.end();it++){
    std::cout << it->first << std::endl;
    for(auto j=0;j<it->second.size();j++){
      std::cout << it->second[j] << " ";
    }
    std::cout <<std::endl;
  }
}
