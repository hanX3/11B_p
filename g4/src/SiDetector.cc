#include "SiDetector.hh"

#include "G4TwoVector.hh"
#include "TString.h"

#include <iostream>
#include <string>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiDetector::SiDetector(G4LogicalVolume* log)
    : exp_hall_log(log)
{
  check_overlaps = true;
}

SiDetector::~SiDetector() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Annular S3 wafer.
void SiDetector::ConstructSiDetector(const std::array<G4double, 2>& si_par,
                                     G4Material* mat)
{
  const TString solid_name = TString::Format("%s_solid", si_name.c_str());
  const auto inner_it = map_si_inner_radius.find(si_name);
  const G4double inner_radius = inner_it != map_si_inner_radius.end()
                                    ? inner_it->second * mm
                                    : 0. * mm;

  auto* solid = new G4Tubs(solid_name.Data(),
                           inner_radius,
                           si_par[0] / 2. * mm,
                           si_par[1] / 2. * mm,
                           0. * deg,
                           360. * deg);

  const TString log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, log_name.Data());

  const auto& colour = map_color_par[si_name];
  auto* vis = new G4VisAttributes(G4Colour(colour[0], colour[1], colour[2], colour[3]));
  vis->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Rectangular W1 wafer.
void SiDetector::ConstructSiBoxDetector(const std::array<G4double, 3>& si_box_par,
                                        G4Material* mat)
{
  const TString solid_name = TString::Format("%s_box_solid", si_name.c_str());
  auto* solid = new G4Box(solid_name.Data(),
                          si_box_par[0] / 2. * mm,
                          si_box_par[1] / 2. * mm,
                          si_box_par[2] / 2. * mm);

  const TString log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, log_name.Data());

  const auto& colour = map_color_par[si_name];
  auto* vis = new G4VisAttributes(G4Colour(colour[0], colour[1], colour[2], colour[3]));
  vis->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Flat trapezoidal wedge DSSSD used by the CAKE-like lampshades.
// Local coordinates:
//   x: transverse/sector direction
//   y: narrow edge (-y) -> wide edge (+y), ring direction
//   z: wafer normal; placement makes local +z face the target
void SiDetector::ConstructSiWedgeDetector(const std::array<G4double, 4>& si_wedge_par,
                                          G4Material* mat)
{
  const G4double inner_half_width = si_wedge_par[0] / 2. * mm;
  const G4double outer_half_width = si_wedge_par[1] / 2. * mm;
  const G4double half_length = si_wedge_par[2] / 2. * mm;
  const G4double half_thickness = si_wedge_par[3] / 2. * mm;

  std::vector<G4TwoVector> polygon = {
      G4TwoVector(-inner_half_width, -half_length),
      G4TwoVector(+inner_half_width, -half_length),
      G4TwoVector(+outer_half_width, +half_length),
      G4TwoVector(-outer_half_width, +half_length)};

  const TString solid_name = TString::Format("%s_wedge_solid", si_name.c_str());
  auto* solid = new G4ExtrudedSolid(solid_name.Data(),
                                    polygon,
                                    half_thickness,
                                    G4TwoVector(),
                                    1.0,
                                    G4TwoVector(),
                                    1.0);

  const TString log_name = TString::Format("%s_log", si_name.c_str());
  si_detector_log = new G4LogicalVolume(solid, mat, log_name.Data());

  const auto& colour = map_color_par[si_name];
  auto* vis = new G4VisAttributes(G4Colour(colour[0], colour[1], colour[2], colour[3]));
  vis->SetForceSolid(true);
  si_detector_log->SetVisAttributes(vis);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiDetector::PlaceSiDetector(G4RotationMatrix* rot, const G4ThreeVector& pos)
{
  const TString physical_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(new G4RotationMatrix(rot->inverse()),
                                      pos,
                                      si_detector_log,
                                      physical_name.Data(),
                                      exp_hall_log,
                                      false,
                                      EncodeDetectorCopyNo(DetectorType::Si,
                                                           0,
                                                           ring_id,
                                                           sector_id),
                                      check_overlaps);
}

void SiDetector::PlaceSiDetector(const G4Transform3D& transform_3d)
{
  const TString physical_name = TString::Format("%s_phy", si_name.c_str());
  si_detector_phy = new G4PVPlacement(transform_3d,
                                      si_detector_log,
                                      physical_name.Data(),
                                      exp_hall_log,
                                      false,
                                      EncodeDetectorCopyNo(DetectorType::Si,
                                                           0,
                                                           ring_id,
                                                           sector_id),
                                      check_overlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Subarray IDs remain stable so existing ROOT analyses can continue to use:
//   1: W1 drum
//   2: backward S3
//   3: forward S3
//   4: forward wedge lampshade  (historical name Si_ForwardCap)
//   5: backward wedge lampshade (historical name Si_BackwardCap)
std::map<G4String, G4int> SiDetector::map_name_to_ring_id = {
    {"Si_Drum", 1},
    {"Si_BackwardAnnular", 2},
    {"Si_ForwardAnnular", 3},
    {"Si_ForwardCap", 4},
    {"Si_BackwardCap", 5}};

std::map<G4int, G4String> SiDetector::map_ring_id_to_name = {
    {1, "Si_Drum"},
    {2, "Si_BackwardAnnular"},
    {3, "Si_ForwardAnnular"},
    {4, "Si_ForwardCap"},
    {5, "Si_BackwardCap"}};

std::map<G4String, G4int> SiDetector::map_name_to_sectors = {
    {"Si_Drum", 12},
    {"Si_BackwardAnnular", 1},
    {"Si_ForwardAnnular", 1},
    {"Si_ForwardCap", 12},
    {"Si_BackwardCap", 12}};

// S3 active outer diameter and thickness [mm].
std::map<G4String, std::array<G4double, 2>> SiDetector::map_si_par = {
    {"Si_BackwardAnnular", {70., 0.5}},
    {"Si_ForwardAnnular", {70., 0.5}}};

// Only the equatorial detector remains a square W1.
std::map<G4String, std::array<G4double, 3>> SiDetector::map_si_box_par = {
    {"Si_Drum", {50., 50., 0.5}}};

// Values are replaced from SiArrayConfig at every geometry construction.
std::map<G4String, std::array<G4double, 4>> SiDetector::map_si_wedge_par = {
    {"Si_ForwardCap", {11.8348065, 50.9890543, 102.5, 0.4}},
    {"Si_BackwardCap", {11.8348065, 50.9890543, 102.5, 0.4}}};

std::map<G4String, G4double> SiDetector::map_si_inner_radius = {
    {"Si_BackwardAnnular", 11.},
    {"Si_ForwardAnnular", 11.}};

// Documentation defaults; SiArray::CalculatePlacement computes live values.
std::map<G4String, std::array<G4double, 3>> SiDetector::map_placement_par = {
    {"Si_Drum", {93.3, 0., TargetZPos / mm}},
    {"Si_BackwardAnnular", {0., 0., TargetZPos / mm - 159.2916073617}},
    {"Si_ForwardAnnular", {0., 0., TargetZPos / mm + 159.2916073617}},
    {"Si_ForwardCap", {62.988, 75.073, TargetZPos / mm + 75.073}},
    {"Si_BackwardCap", {62.988, 75.073, TargetZPos / mm - 75.073}}};

std::map<G4String, std::array<G4double, 4>> SiDetector::map_color_par = {
    {"Si_Drum", {0., 0., 1., 0.8}},
    {"Si_ForwardCap", {0.15, 0.75, 0.95, 0.8}},
    {"Si_BackwardCap", {0.65, 0.2, 0.95, 0.8}},
    {"Si_BackwardAnnular", {0., 0.45, 1., 0.8}},
    {"Si_ForwardAnnular", {0., 1., 0.45, 0.8}}};
