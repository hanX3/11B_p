#include "SiDetector.hh"

#include <string>
#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiDetector::SiDetector(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  check_overlaps = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiDetector::~SiDetector() {}

// Si
void SiDetector::ConstructSiDetector(const std::array<G4double, 2>& si_par, G4Material* mat)
{
  std::cout << "start const Si Detector." << std::endl;

  TString si_solid_name = TString::Format("%s_solid", si_name.c_str());
  const auto inner_it = SiDetector::map_si_inner_radius.find(si_name);
  const G4double inner_radius = inner_it != SiDetector::map_si_inner_radius.end() ? inner_it->second * mm : 0. * mm;
  G4Tubs* solid = new G4Tubs(si_solid_name.Data(), inner_radius, si_par[0] / 2. * mm, si_par[1] / 2. * mm, 0. * deg, 360. * deg);

  TString si_log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, si_log_name.Data());

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(SiDetector::map_color_par[si_name][0], SiDetector::map_color_par[si_name][1], SiDetector::map_color_par[si_name][2], SiDetector::map_color_par[si_name][3]));
  vis_att->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis_att);
}

// Si rectangular barrel module
void SiDetector::ConstructSiBoxDetector(const std::array<G4double, 3>& si_box_par, G4Material* mat)
{
  std::cout << "start const Si box Detector." << std::endl;

  TString si_solid_name = TString::Format("%s_box_solid", si_name.c_str());
  G4Box* solid = new G4Box(si_solid_name.Data(), si_box_par[0] / 2. * mm, si_box_par[1] / 2. * mm, si_box_par[2] / 2. * mm);

  TString si_log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, si_log_name.Data());

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(SiDetector::map_color_par[si_name][0], SiDetector::map_color_par[si_name][1], SiDetector::map_color_par[si_name][2], SiDetector::map_color_par[si_name][3]));
  vis_att->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiDetector::PlaceSiDetector(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString si_phy_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, si_detector_log, si_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::Si, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiDetector::PlaceSiDetector(const G4Transform3D& transfrom_3d)
{
  TString si_phy_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(transfrom_3d, si_detector_log, si_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::Si, 0, ring_id, sector_id), check_overlaps);
}

// Al shell
void SiDetector::ConstructAlShell(const std::array<G4double, 6>& al_par, G4Material* mat)
{
  std::cout << "start const Al Shell Detector." << std::endl;

  G4double zz_z_plane[4] = {0., al_par[4] * mm, al_par[4] * mm, al_par[5] * mm};
  G4double zz_r_inner_plane[4] = {al_par[0] * mm, al_par[0] * mm, al_par[2] * mm, al_par[2] * mm};
  G4double zz_r_outer_plane[4] = {al_par[1] * mm, al_par[1] * mm, al_par[3] * mm, al_par[3] * mm};

  G4Polycone* solid = new G4Polycone("solid", 0. * deg, 360. * deg, 4, zz_z_plane, zz_r_inner_plane, zz_r_outer_plane);

  TString al_shell_log_name = TString::Format("%s_al_shell_log", si_name.c_str());
  al_shell_log = new G4LogicalVolume(solid, mat, al_shell_log_name.Data());

  TString al_shell_reg_name = TString::Format("%s_%d_SiAlShell", si_name.c_str(), sector_id);
  G4Region* al_shell_reg = new G4Region(al_shell_reg_name.Data());
  al_shell_reg->AddRootLogicalVolume(al_shell_log);

  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(SiDetector::map_color_al_shell_par[si_name][0], SiDetector::map_color_al_shell_par[si_name][1], SiDetector::map_color_al_shell_par[si_name][2], SiDetector::map_color_al_shell_par[si_name][3]));
  vis_att->SetForceSolid(true);
  al_shell_log->SetVisAttributes(vis_att);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiDetector::PlaceAlShell(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  TString al_shell_phy_name = TString::Format("%s_al_shell_phy", si_name.c_str());
  al_shell_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()), pos, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::Si, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiDetector::PlaceAlShell(const G4Transform3D& transfrom_3d)
{
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
  al_shell_phy = new G4PVPlacement(transform_3d_new, al_shell_log, al_shell_phy_name.Data(), exp_hall_log, false, EncodeDetectorCopyNo(DetectorType::Si, 0, ring_id, sector_id), check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// All-commercial "connected drum" baseline (W1 + S3 catalogue parts).
//   ring_id 1: Si_Drum            - 12 W1 faces, dodecagonal prism, theta ~75-105
//   ring_id 2: Si_BackwardAnnular - S3 annular DSSD at -36 mm, theta ~136-163
//   ring_id 3: Si_ForwardAnnular  - S3 annular DSSD at +36 mm, theta ~17-44
//   ring_id 4: Si_ForwardCap      - 12 W1 fish-scale plates hinged on the drum
//                                   forward edge, folded inward, theta ~44-75
//   ring_id 5: Si_BackwardCap     - mirror of ring 4, theta ~105-136
// Even cap module ids form the inner fish-scale sub-ring (hinged on the drum
// edge); odd ids form the outer sub-ring (capStagger further out, rotated half
// a face period) covering the V-gaps.
std::map<G4String, G4int> SiDetector::map_name_to_ring_id = {{"Si_Drum", 1}, {"Si_BackwardAnnular", 2}, {"Si_ForwardAnnular", 3}, {"Si_ForwardCap", 4}, {"Si_BackwardCap", 5}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4int, G4String> SiDetector::map_ring_id_to_name = {{1, "Si_Drum"}, {2, "Si_BackwardAnnular"}, {3, "Si_ForwardAnnular"}, {4, "Si_ForwardCap"}, {5, "Si_BackwardCap"}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Drum/cap module counts are documented defaults; the value actually used at
// build time is SiArrayConfig::GetDrumModules() so it can be scanned.
std::map<G4String, G4int> SiDetector::map_name_to_sectors = {{"Si_Drum", 12}, {"Si_BackwardAnnular", 1}, {"Si_ForwardAnnular", 1}, {"Si_ForwardCap", 12}, {"Si_BackwardCap", 12}};

// 0: diameter (S3 active outer diameter 70 mm)
// 1: thickness
std::map<G4String, std::array<G4double, 2>> SiDetector::map_si_par = {{"Si_BackwardAnnular", {70., 0.5}}, {"Si_ForwardAnnular", {70., 0.5}}};

// W1-type quadratic DSSD: 50 x 50 mm active area, 0.5 mm thick.
// 0: width  (azimuthal direction, local x)
// 1: length (polar direction,   local y)
// 2: thickness
std::map<G4String, std::array<G4double, 3>> SiDetector::map_si_box_par = {{"Si_Drum", {50., 50., 0.5}}, {"Si_ForwardCap", {50., 50., 0.5}}, {"Si_BackwardCap", {50., 50., 0.5}}};

// S3 active inner radius 11 mm.
std::map<G4String, G4double> SiDetector::map_si_inner_radius = {{"Si_BackwardAnnular", 11.}, {"Si_ForwardAnnular", 11.}};

// 0: x
// 1: y
// 2: z
// Documented defaults only.  The actual placement of every sub-array is
// computed inside SiArray::CalculatePlacement from SiArrayConfig (drum
// inscribed radius derived from the module count; cap fold angle / stagger;
// annular distances), so everything can be scanned from a macro.
std::map<G4String, std::array<G4double, 3>> SiDetector::map_placement_par = {{"Si_Drum", {93.3, 0., TargetZPos / mm}}, {"Si_BackwardAnnular", {0., 0., TargetZPos / mm - 36.}}, {"Si_ForwardAnnular", {0., 0., TargetZPos / mm + 36.}}, {"Si_ForwardCap", {93.3, 45., TargetZPos / mm + 25.}}, {"Si_BackwardCap", {93.3, 45., TargetZPos / mm - 25.}}};

// 0: r1_inner
// 1: r1_outer
// 2: r2_inner
// 3: r2_outer
// 4: h_r1
// 5: h_r2
// Resized for the S3 active outer radius (35 mm): collar 36-37 mm plus a
// backing plate down to r = 9 mm (physical hole slightly inside the 11 mm
// active inner radius).
std::map<G4String, std::array<G4double, 6>> SiDetector::map_al_par = {{"Si_BackwardAnnular", {36., 37., 9., 37., 0.5, 1.5}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> SiDetector::map_color_par = {{"Si_Drum", {0., 0., 1., 0.8}}, {"Si_ForwardCap", {0.3, 0., 1., 0.8}}, {"Si_BackwardCap", {0.5, 0., 1., 0.8}}, {"Si_BackwardAnnular", {0., 0.45, 1., 0.8}}, {"Si_ForwardAnnular", {0., 1., 0.45, 0.8}}};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::map<G4String, std::array<G4double, 4>> SiDetector::map_color_al_shell_par = {{"Si_BackwardAnnular", {0.8, 0., 1., 0.3}}};
