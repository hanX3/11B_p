#include "LaBr3Detector.hh"

#include <string>
#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Detector::LaBr3Detector(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  check_overlaps = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Detector::~LaBr3Detector() {}

// LaBr3
void LaBr3Detector::ConstructLaBr3Detector(const std::array<G4double, 2>& labr3_par, G4Material* mat)
{
  std::cout << "start const LaBr3 Detector." << std::endl;

  G4Tubs* labr3_detector_solid = new G4Tubs("labr3_detector_solid", 0. * mm, labr3_par[0] / 2. * mm, labr3_par[1] / 2. * mm, 0. * deg, 360. * deg);

  TString labr3_log_name = TString::Format("%s_log", labr3_name.c_str());
  labr3_detector_log = new G4LogicalVolume(labr3_detector_solid, mat, labr3_log_name.Data());

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(LaBr3Detector::map_color_par[labr3_name][0], LaBr3Detector::map_color_par[labr3_name][1], LaBr3Detector::map_color_par[labr3_name][2], LaBr3Detector::map_color_par[labr3_name][3]));
  vis_att->SetForceSolid(true);
  labr3_detector_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Detector::PlaceLaBr3Detector(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString labr3_phy_name = TString::Format("%s_phy", labr3_name.c_str());
  labr3_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, labr3_detector_log, labr3_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::LaBr3, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Detector::PlaceLaBr3Detector(const G4Transform3D& transfrom_3d)
{
  TString labr3_phy_name = TString::Format("%s_phy", labr3_name.c_str());
  labr3_detector_phy = new G4PVPlacement(transfrom_3d, labr3_detector_log, labr3_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::LaBr3, 0, ring_id, sector_id), check_overlaps);
}

// Al shell
void LaBr3Detector::ConstructAlShell(const std::array<G4double, 3>& al_par, G4Material* mat)
{
  std::cout << "start const Al Shell Detector." << std::endl;

  G4Tubs* solid = new G4Tubs("al_shell_solid", al_par[0] / 2. * mm, al_par[1] / 2. * mm, al_par[2] / 2. * mm, 0. * deg, 360. * deg);

  TString al_shell_log_name = TString::Format("%s_al_shell_log", labr3_name.c_str());
  al_shell_log = new G4LogicalVolume(solid, mat, al_shell_log_name.Data());

  TString al_shell_reg_name = TString::Format("%s_%d_LaBr3AlShell", labr3_name.c_str(), sector_id);
  G4Region* al_shell_reg = new G4Region(al_shell_reg_name.Data());
  al_shell_reg->AddRootLogicalVolume(al_shell_log);

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(LaBr3Detector::map_color_al_shell_par[labr3_name][0], LaBr3Detector::map_color_al_shell_par[labr3_name][1], LaBr3Detector::map_color_al_shell_par[labr3_name][2], LaBr3Detector::map_color_al_shell_par[labr3_name][3]));
  vis_att->SetForceSolid(true);
  al_shell_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Detector::PlaceAlShell(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", labr3_name.c_str());
  al_shell_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::LaBr3, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Detector::PlaceAlShell(const G4Transform3D& transfrom_3d)
{
  G4cout << "----> start place Al shell" << G4endl;

  G4RotationMatrix rot = transfrom_3d.getRotation();
  rot.print(G4cout);
  G4ThreeVector pos = transfrom_3d.getTranslation();

  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  G4cout << "x = " << x << " y = " << y << " z = " << z << G4endl;

  G4ThreeVector pos_new = G4ThreeVector(x, y, z);
  G4Transform3D transform_3d_new(rot, pos_new);

  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", labr3_name.c_str());
  al_shell_phy = new G4PVPlacement(transform_3d_new, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::LaBr3, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, G4int> LaBr3Detector::map_name_to_ring_id = {{"LaBr3_NearSide", 1}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4int, G4String> LaBr3Detector::map_ring_id_to_name = {{1, "LaBr3_NearSide"}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, G4int> LaBr3Detector::map_name_to_sectors = {{"LaBr3_NearSide", 4}};

// 0: diameter
// 1: length
std::map<G4String, std::array<G4double, 2>> LaBr3Detector::map_labr3_par = {{"LaBr3_NearSide", {50.8, 76.2}}};

// 0: x
// 1: y
// 2: z
std::map<G4String, std::array<G4double, 3>> LaBr3Detector::map_placement_par = {{"LaBr3_NearSide", {225.6, 0., 158.5}}};

// 0: diameter inner
// 1: diameter outer
// 1: length
std::map<G4String, std::array<G4double, 3>> LaBr3Detector::map_al_par = {{"LaBr3_NearSide", {51.8, 55.8, 78.2}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> LaBr3Detector::map_color_par = {{"LaBr3_NearSide", {0.4, 0., 0.8, 0.3}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> LaBr3Detector::map_color_al_shell_par = {{"LaBr3_NearSide", {0.2, 0., 1., 0.3}}};
