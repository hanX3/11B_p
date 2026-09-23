#include "HPGeDetector.hh"

#include <string>
#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HPGeDetector::HPGeDetector(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  check_overlaps = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HPGeDetector::~HPGeDetector() {}

// HPGe
void HPGeDetector::ConstructHPGeDetector(const std::array<G4double, 4>& hpge_par, G4Material* mat)
{
  std::cout << "start const HPGe Detector." << std::endl;

  G4Tubs* hpge_crystal_solid = new G4Tubs("hpge_crystal_solid", 0. * mm, hpge_par[0] / 2. * mm, hpge_par[1] / 2. * mm, 0. * deg, 360. * deg);
  G4Tubs* hpge_ln2_solid = new G4Tubs("hpge_hole_solid", 0. * mm, hpge_par[2] / 2. * mm, hpge_par[3] / 2. * mm, 0. * deg, 360. * deg);
  G4double zz = hpge_par[1] * mm - hpge_par[3] * mm;
  G4SubtractionSolid* hpge_detector_solid = new G4SubtractionSolid("hpge_detector_solid", hpge_crystal_solid, hpge_ln2_solid, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, zz)));

  TString hpge_log_name = TString::Format("%s_log", hpge_name.c_str());
  hpge_detector_log = new G4LogicalVolume(hpge_detector_solid, mat, hpge_log_name.Data());

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(HPGeDetector::map_color_par[hpge_name][0], HPGeDetector::map_color_par[hpge_name][1], HPGeDetector::map_color_par[hpge_name][2], HPGeDetector::map_color_par[hpge_name][3]));
  vis_att->SetForceSolid(true);
  hpge_detector_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeDetector::PlaceHPGeDetector(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString hpge_phy_name = TString::Format("%s_phy", hpge_name.c_str());
  hpge_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, hpge_detector_log, hpge_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::HPGe, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeDetector::PlaceHPGeDetector(const G4Transform3D& transfrom_3d)
{
  TString hpge_phy_name = TString::Format("%s_phy", hpge_name.c_str());
  hpge_detector_phy = new G4PVPlacement(transfrom_3d, hpge_detector_log, hpge_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::HPGe, 0, ring_id, sector_id), check_overlaps);
}

// Al shell
void HPGeDetector::ConstructAlShell(const std::array<G4double, 3>& al_par, G4Material* mat)
{
  std::cout << "start const Al Shell Detector." << std::endl;

  G4Tubs* solid = new G4Tubs("al_shell_solid", al_par[0] / 2. * mm, al_par[1] / 2. * mm, al_par[2] / 2. * mm, 0. * deg, 360. * deg);

  TString al_shell_log_name = TString::Format("%s_al_shell_log", hpge_name.c_str());
  al_shell_log = new G4LogicalVolume(solid, mat, al_shell_log_name.Data());

  TString al_shell_reg_name = TString::Format("%s_%d_HPGeAlShell", hpge_name.c_str(), sector_id);
  G4Region* al_shell_reg = new G4Region(al_shell_reg_name.Data());
  al_shell_reg->AddRootLogicalVolume(al_shell_log);

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(HPGeDetector::map_color_al_shell_par[hpge_name][0], HPGeDetector::map_color_al_shell_par[hpge_name][1], HPGeDetector::map_color_al_shell_par[hpge_name][2], HPGeDetector::map_color_al_shell_par[hpge_name][3]));
  vis_att->SetForceSolid(true);
  al_shell_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeDetector::PlaceAlShell(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", hpge_name.c_str());
  al_shell_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::HPGe, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeDetector::PlaceAlShell(const G4Transform3D& transfrom_3d)
{
  G4RotationMatrix rot = transfrom_3d.getRotation();
  rot.print(G4cout);
  G4ThreeVector pos = rot.inverse() * transfrom_3d.getTranslation();

  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  G4cout << "x = " << x << " y = " << y << " z = " << z << G4endl;

  G4ThreeVector pos_new = rot * G4ThreeVector(x, y, z);
  G4Transform3D transform_3d_new(rot, pos_new);

  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", hpge_name.c_str());
  al_shell_phy = new G4PVPlacement(transform_3d_new, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::HPGe, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
// Transverse center radius for the 45/135-degree HPGe rings.
// 240 mm preserves target alignment while moving the HPGe housings outside
// the chamber flange envelope.
constexpr G4double HPGeRingTransverseRadiusMm = 240.0;
}

std::map<G4String, G4int> HPGeDetector::map_name_to_ring_id = {{"HPGe_Forward45", 2}, {"HPGe_Backward135", 3}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4int, G4String> HPGeDetector::map_ring_id_to_name = {{2, "HPGe_Forward45"}, {3, "HPGe_Backward135"}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, G4int> HPGeDetector::map_name_to_sectors = {{"HPGe_Forward45", 4}, {"HPGe_Backward135", 4}};

// 0: diameter
// 1: length
// 2: LN2 hole diameter
// 3: LN2 length
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_hpge_par = {
    {"HPGe_Forward45", {93.4, 54.5, 11.1, 47.95}},
    {"HPGe_Backward135", {93.4, 54.5, 11.1, 47.95}}};

// 0: transverse x reference in mm
// 1: transverse y reference in mm
// 2: absolute z position in mm
//
// The four sectors are placed at azimuths 0, 180, 90, and 270 deg by HPGeArray::SidePortPosition().
// The detector axes are aimed at the target position. With TargetZPos = 170 mm and
// transverse radius = 240 mm, these positions correspond to polar angles of 45 deg
// and 135 deg relative to the beam axis at the target.
std::map<G4String, std::array<G4double, 3>> HPGeDetector::map_placement_par = {
    {"HPGe_Forward45", {HPGeRingTransverseRadiusMm, 0., TargetZPos / mm + HPGeRingTransverseRadiusMm}},
    {"HPGe_Backward135", {HPGeRingTransverseRadiusMm, 0., TargetZPos / mm - HPGeRingTransverseRadiusMm}}};

// 0: diameter inner
// 1: diameter outer
// 2: length
std::map<G4String, std::array<G4double, 3>> HPGeDetector::map_al_par = {
    {"HPGe_Forward45", {94.4, 98.4, 56.5}},
    {"HPGe_Backward135", {94.4, 98.4, 56.5}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_color_par = {
    {"HPGe_Forward45", {0., 1., 1., 0.3}},
    {"HPGe_Backward135", {0., 0.8, 1., 0.3}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> HPGeDetector::map_color_al_shell_par = {
    {"HPGe_Forward45", {0.2, 0.6, 1., 0.3}},
    {"HPGe_Backward135", {0.2, 0.4, 1., 0.3}}};
