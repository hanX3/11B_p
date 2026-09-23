#include "SiArray.hh"
#include "SiDetector.hh"

#include <cmath>
#include "TString.h"

//
SiArray::SiArray(G4LogicalVolume *log)
: exp_hall_log(log)
{
  si_numbers = 0;

  for(auto it=SiDetector::map_name_to_sectors.begin();it!=SiDetector::map_name_to_sectors.end();it++){
    si_numbers += it->second;
  }
  
  //
  G4NistManager *nist_manager = G4NistManager::Instance();

  si_mat = nist_manager->FindOrBuildMaterial("G4_Si");
  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");
  
  PrintDetectorDimensionInfo();
}

//
SiArray::~SiArray()
{

}

//
void SiArray::Construct()
{
  std::vector<SiDetector*>::iterator it = v_si_detector.begin();
  // clear all elements from the array
  for(;it!=v_si_detector.end();it++) delete *it;
  v_si_detector.clear();

  std::map<G4int, G4int> map_i2ring;
  std::map<G4int, G4int> map_i2sector;
  std::map<G4int, G4String> map_i2name;
  int ii = 0;
  for(auto it=SiDetector::map_name_to_ring_id.begin();it!=SiDetector::map_name_to_ring_id.end();it++){
    for(int j=0;j<SiDetector::map_name_to_sectors[it->first];j++){
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

  for(int i=0;i<si_numbers;i++){
    v_si_detector.push_back(new SiDetector(exp_hall_log));
    v_si_detector[i]->SetName(map_i2name[i]);
    v_si_detector[i]->SetRingId(map_i2ring[i]);
    v_si_detector[i]->SetSectorId(map_i2sector[i]);
  }

  for(it=v_si_detector.begin();it!=v_si_detector.end();it++){
    (*it)->ConstructSiDetector(SiDetector::map_si_par[(*it)->GetName()], si_mat);
    (*it)->PlaceSiDetector(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));

    (*it)->ConstructAlShell(SiDetector::map_al_par[(*it)->GetName()], al_mat);
    (*it)->PlaceAlShell(CalculatePlacement((*it)->GetName(), (*it)->GetSectorId()));
  }
}

//
void SiArray::MakeSensitive(SiSD *si_sd)
{
  std::vector<SiDetector*>::iterator it = v_si_detector.begin();
  for(;it!=v_si_detector.end();it++){
     (*it)->GetLog()->SetSensitiveDetector(si_sd);
  }
}

//
G4Transform3D SiArray::CalculatePlacement(G4String name, G4int sector_id)
{
  G4RotationMatrix *rot_matrix = new G4RotationMatrix();
  G4double rot_x_angle = 0. *deg;
  G4double rot_y_angle = 180. *deg;
  G4double rot_z_angle = 0. *deg;
  G4cout << "x angle " << rot_x_angle << " y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  rot_matrix->rotateX(rot_x_angle);
  rot_matrix->rotateY(rot_y_angle);
  rot_matrix->rotateZ(rot_z_angle);

  rot_matrix->print(G4cout);

  G4double x = SiDetector::map_placement_par[name][0] *mm;
  G4double y = SiDetector::map_placement_par[name][1] *mm;
  G4double z = SiDetector::map_placement_par[name][2] *mm;
  G4ThreeVector pos = G4ThreeVector(x, y, z);
  G4cout << "pos x " << pos.x() << "pos y " << pos.y() << "pos z " << pos.z() << G4endl;

  G4Transform3D transform(*rot_matrix, pos);

  return transform;
}

//
void SiArray::PrintDetectorDimensionInfo()
{
  for(auto it=SiDetector::map_si_par.begin();it!=SiDetector::map_si_par.end();it++){
    std::cout << it->first << std::endl;
    for(auto j=0;j<it->second.size();j++){
      std::cout << it->second[j] << " ";
    }
    std::cout <<std::endl;
  }
}
