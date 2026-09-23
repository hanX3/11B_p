#include "SiDetector.hh"

#include <string>
#include "TString.h"

//
SiDetector::SiDetector(G4LogicalVolume *log)
: exp_hall_log(log)
{
  check_overlaps = true;
}

//
SiDetector::~SiDetector()
{

}

// Si
void SiDetector::ConstructSiDetector(const std::array<G4double, 2> &si_par, G4Material *mat)
{
  std::cout << "start const Si Detector." << std::endl;

  TString si_solid_name = TString::Format("%s_solid", si_name.c_str());
  // std::cout << si_solid_name << std::endl;
  
  G4Tubs *solid = new G4Tubs("solid", 0.*mm, SiDetector::map_si_par[si_name][0]/2. *mm, SiDetector::map_si_par[si_name][1]/2. *mm, 0. *deg, 360. *deg);

  TString si_log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, si_log_name.Data());
  
  // color
  G4VisAttributes *vis_att = new G4VisAttributes(G4Colour(SiDetector::map_color_par[si_name][0], SiDetector::map_color_par[si_name][1], SiDetector::map_color_par[si_name][2], SiDetector::map_color_par[si_name][3]));
  vis_att->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis_att);
}

//
void SiDetector::PlaceSiDetector(G4RotationMatrix *rot, const G4ThreeVector &pos)
{
  TString si_phy_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, si_detector_log, si_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

//
void SiDetector::PlaceSiDetector(const G4Transform3D &transfrom_3d)
{
  TString si_phy_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(transfrom_3d, si_detector_log, si_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

// Al shell
void SiDetector::ConstructAlShell(const std::array<G4double, 6> &al_par, G4Material *mat)
{
  std::cout << "start const Al Shell Detector." << std::endl;

  G4double zz_z_plane[4] = {0., SiDetector::map_al_par[si_name][4] *mm, SiDetector::map_al_par[si_name][4] *mm, SiDetector::map_al_par[si_name][5] *mm};
  G4double zz_r_inner_plane[4] = {SiDetector::map_al_par[si_name][0] *mm, SiDetector::map_al_par[si_name][0] *mm, SiDetector::map_al_par[si_name][2] *mm, SiDetector::map_al_par[si_name][2] *mm};
  G4double zz_r_outer_plane[4] = {SiDetector::map_al_par[si_name][1] *mm, SiDetector::map_al_par[si_name][1] *mm, SiDetector::map_al_par[si_name][3] *mm, SiDetector::map_al_par[si_name][3] *mm};
  
  G4Polycone *solid = new G4Polycone("solid", 0. *deg, 360. *deg, 4, zz_z_plane, zz_r_inner_plane, zz_r_outer_plane);

  TString al_shell_log_name = TString::Format("%s_al_shell_log", si_name.c_str());
  al_shell_log = new G4LogicalVolume(solid, mat, al_shell_log_name.Data());
  
  // G4Region for Cuts
  G4Region *al_shell_reg = new G4Region("SiAlShell");
  al_shell_reg->AddRootLogicalVolume(al_shell_log);

  // color
  G4VisAttributes *vis_att = new G4VisAttributes(G4Colour(SiDetector::map_color_al_shell_par[si_name][0], 
                                                          SiDetector::map_color_al_shell_par[si_name][1], 
                                                          SiDetector::map_color_al_shell_par[si_name][2], 
                                                          SiDetector::map_color_al_shell_par[si_name][3]));
  vis_att->SetForceSolid(true);
  al_shell_log->SetVisAttributes(vis_att);
}

//
void SiDetector::PlaceAlShell(G4RotationMatrix *rot, const G4ThreeVector &pos)
{
  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", si_name.c_str());
  al_shell_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

//
void SiDetector::PlaceAlShell(const G4Transform3D &transfrom_3d)
{
  // adjust transfrom_3d
  G4RotationMatrix rot = transfrom_3d.getRotation();
  rot.print(G4cout);

  G4ThreeVector pos = transfrom_3d.getTranslation();

  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  G4cout << "x = " << x << " y = " << y << " z = " << z << G4endl;

  G4ThreeVector pos_new = G4ThreeVector(x, y, z);
  G4Transform3D transform_3d_new(rot, pos_new);

  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(transform_3d_new, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}


//
std::map<G4String, G4int> SiDetector::map_name_to_ring_id = {
  {"Si_01", 1}
};

//
std::map<G4int, G4String> SiDetector::map_ring_id_to_name = {
  {1, "Si_01"}
};

//
std::map<G4String, G4int> SiDetector::map_name_to_sectors = {
  {"Si_01", 1}
};

// 0: diameter
// 1: length
std::map<G4String, std::array<G4double, 2>> SiDetector::map_si_par = {
  {"Si_01", {40., 0.3}}
};

// 0: x
// 1: y
// 2: z
std::map<G4String, std::array<G4double, 3>> SiDetector::map_placement_par = {
  {"Si_01", {105.2, 0, -167.}}
};

// 0: r1_inner
// 1: r1_outer
// 2: r2_inner
// 3: r2_outer
// 4: h_r1
// 5: h_r2
std::map<G4String, std::array<G4double, 6>> SiDetector::map_al_par = {
  {"Si_01", {22., 25, 0., 25, 0.5, 11.5}}
};

//
std::map<G4String, std::array<G4double, 4>> SiDetector::map_color_par = {
  {"Si_01", {0., 0., 1., 0.8}} // blue
};

//
std::map<G4String, std::array<G4double, 4>> SiDetector::map_color_al_shell_par = {
  {"Si_01", {0.8, 0., 1., 0.3}}
};
