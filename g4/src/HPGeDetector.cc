#include "HPGeDetector.hh"

#include <string>
#include "TString.h"

//
HPGeDetector::HPGeDetector(G4LogicalVolume *log)
: exp_hall_log(log)
{
  check_overlaps = true;
}

//
HPGeDetector::~HPGeDetector()
{

}

// HPGe
void HPGeDetector::ConstructHPGeDetector(const std::array<G4double, 4> &hpge_par, G4Material *mat)
{
  std::cout << "start const HPGe Detector." << std::endl;
  
  G4Tubs *hpge_crystal_solid = new G4Tubs("hpge_crystal_solid", 0.*mm, 
                                   HPGeDetector::map_hpge_par[hpge_name][0]/2. *mm, 
                                   HPGeDetector::map_hpge_par[hpge_name][1]/2. *mm,
                                   0. *deg,
                                   360. *deg);

  G4Tubs *hpge_ln2_solid = new G4Tubs("hpge_hole_solid", 0.*mm, 
                               HPGeDetector::map_hpge_par[hpge_name][2]/2. *mm, 
                               HPGeDetector::map_hpge_par[hpge_name][3]/2. *mm,
                               0. *deg,
                               360. *deg);

  G4double zz = HPGeDetector::map_hpge_par[hpge_name][1]*mm-HPGeDetector::map_hpge_par[hpge_name][3]*mm;
  G4SubtractionSolid *hpge_detector_solid = new G4SubtractionSolid("hpge_detector_solid",
                                                hpge_crystal_solid, 
                                                hpge_ln2_solid, 
                                                G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, zz)));

  TString hpge_log_name = TString::Format("%s_log", hpge_name.c_str());
  hpge_detector_log = new G4LogicalVolume(hpge_detector_solid, mat, hpge_log_name.Data());
  
  // color
  G4VisAttributes *vis_att = new G4VisAttributes(G4Colour(HPGeDetector::map_color_par[hpge_name][0], 
                                 HPGeDetector::map_color_par[hpge_name][1], 
                                 HPGeDetector::map_color_par[hpge_name][2], 
                                 HPGeDetector::map_color_par[hpge_name][3]));
  vis_att->SetForceSolid(true);
  hpge_detector_log->SetVisAttributes(vis_att);
}

//
void HPGeDetector::PlaceHPGeDetector(G4RotationMatrix *rot, const G4ThreeVector &pos)
{
  TString hpge_phy_name = TString::Format("%s_phy", hpge_name.c_str());
  hpge_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, hpge_detector_log, hpge_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

//
void HPGeDetector::PlaceHPGeDetector(const G4Transform3D &transfrom_3d)
{
  TString hpge_phy_name = TString::Format("%s_phy", hpge_name.c_str());
  hpge_detector_phy = new G4PVPlacement(transfrom_3d, hpge_detector_log, hpge_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

// Al shell
void HPGeDetector::ConstructAlShell(const std::array<G4double, 3> &al_par, G4Material *mat)
{
  std::cout << "start const Al Shell Detector." << std::endl;

  G4Tubs *solid = new G4Tubs("al_shell_solid",
                             HPGeDetector::map_al_par[hpge_name][0]/2. *mm, 
                             HPGeDetector::map_al_par[hpge_name][1]/2. *mm,
                             HPGeDetector::map_al_par[hpge_name][2]/2. *mm,
                             0. *deg,
                             360. *deg);

  TString al_shell_log_name = TString::Format("%s_al_shell_log", hpge_name.c_str());
  al_shell_log = new G4LogicalVolume(solid, mat, al_shell_log_name.Data());

  // G4Region for Cuts
  G4Region *al_shell_reg = new G4Region("HPGeAlShell");
  al_shell_reg->AddRootLogicalVolume(al_shell_log);
  
  // color
  G4VisAttributes *vis_att = new G4VisAttributes(G4Colour(HPGeDetector::map_color_al_shell_par[hpge_name][0], 
                                                          HPGeDetector::map_color_al_shell_par[hpge_name][1], 
                                                          HPGeDetector::map_color_al_shell_par[hpge_name][2], 
                                                          HPGeDetector::map_color_al_shell_par[hpge_name][3]));
  vis_att->SetForceSolid(true);
  al_shell_log->SetVisAttributes(vis_att);
}

//
void HPGeDetector::PlaceAlShell(G4RotationMatrix *rot, const G4ThreeVector &pos)
{
  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", hpge_name.c_str());
  al_shell_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

//
void HPGeDetector::PlaceAlShell(const G4Transform3D &transfrom_3d)
{
  // adjust transfrom_3d
  G4RotationMatrix rot = transfrom_3d.getRotation();
  rot.print(G4cout);
  G4ThreeVector pos = rot.inverse()*transfrom_3d.getTranslation();

  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  G4cout << "x = " << x << " y = " << y << " z = " << z << G4endl;

  G4ThreeVector pos_new = rot*G4ThreeVector(x, y, z);
  G4Transform3D transform_3d_new(rot, pos_new);

  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", hpge_name.c_str());
  hpge_detector_phy = new G4PVPlacement(transform_3d_new, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, sector_id, check_overlaps);
}

//
std::map<G4String, G4int> HPGeDetector::map_name_to_ring_id = {
  {"HPGe_01", 1}
};

//
std::map<G4int, G4String> HPGeDetector::map_ring_id_to_name = {
  {1, "HPGe_01"}
};

//
std::map<G4String, G4int> HPGeDetector::map_name_to_sectors = {
  {"HPGe_01", 1}
};

// 0: diameter
// 1: length
// 2: LN2 hole diameter
// 3: LN2 length
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_hpge_par = {
  {"HPGe_01", {93.4, 54.5, 11.1, 47.95}}
};

// 0: x
// 1: y
// 2: z
std::map<G4String, std::array<G4double, 3>> HPGeDetector::map_placement_par = {
  {"HPGe_01", {-214.75, 0., 158.5}}
};

// 0: diameter inner
// 1: diameter outer
// 1: length
std::map<G4String, std::array<G4double, 3>> HPGeDetector::map_al_par = {
  {"HPGe_01", {94.4, 98.4, 56.5}}
};

//
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_color_par = {
  {"HPGe_01", {0., 1., 1., 0.3}}
};

//
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_color_al_shell_par = {
  {"HPGe_01", {0.2, 0.6, 1., 0.3}}
};
