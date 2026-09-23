#include "DetectorConstruction.hh"

#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"
#include "G4GenericMessenger.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "VirtualSphereConfig.hh"
#include "VirtualSphereSD.hh"

#include "G4SystemOfUnits.hh"

#include <string>
#include "TString.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
{
  si_array = nullptr;
  virtual_sphere_log = nullptr;
  hpge_array = nullptr;
  labr3_array = nullptr;

  air_mat = nullptr;
  vaccum_mat = nullptr;
  al_mat = nullptr;
  si_mat = nullptr;
  stainless_steel_mat = nullptr;
  enriched_11b_mat = nullptr;
  natured_11b_high_density_mat = nullptr;
  natured_11b_low_density_mat = nullptr;
  borated_pe_mat = nullptr;
  hb_mat = nullptr;

  DefineMaterials();

  check_overlaps = true;

  //
  SetTargetThickness(TargetThickness);
  SetTargetMaterial(TargetMaterial);

  SetTargetBackingFlag(TargetBackingFlag);

  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  delete step_limit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nist_manager = G4NistManager::Instance();

  air_mat = nist_manager->FindOrBuildMaterial("G4_AIR");
  vaccum_mat = nist_manager->FindOrBuildMaterial("G4_Galactic");
  al_mat = nist_manager->FindOrBuildMaterial("G4_Al");
  si_mat = nist_manager->FindOrBuildMaterial("G4_Si");
  stainless_steel_mat = nist_manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  // Enriched 11B
  G4Isotope* b10 = new G4Isotope("B10", 5, 10, 10 * g / mole);
  G4Isotope* b11 = new G4Isotope("B11", 5, 11, 11 * g / mole);
  G4Element* enriched_b_el = new G4Element("Boron", "B", 2);
  enriched_b_el->AddIsotope(b10, 0.01);
  enriched_b_el->AddIsotope(b11, 0.99);
  enriched_11b_mat = new G4Material("Enriched_B11", 2.38 * g / cm3, 1);
  enriched_11b_mat->AddElement(enriched_b_el, 1);

  G4Element* natured_b_el = new G4Element("Boron", "B", 2);
  natured_b_el->AddIsotope(b10, 0.2);
  natured_b_el->AddIsotope(b11, 0.8);

  // Natured 11B high density
  natured_11b_high_density_mat = new G4Material("Natured_11B_high_density", 2.31 * g / cm3, 1);
  natured_11b_high_density_mat->AddElement(natured_b_el, 1);

  // Natured 11B low density
  natured_11b_low_density_mat = new G4Material("Natured_11B_low_density", 1.4 * g / cm3, 1);
  natured_11b_low_density_mat->AddElement(natured_b_el, 1);

  // Borated polyethylene
  borated_pe_mat = new G4Material("Borated_polyethylene", 1.13 * g / cm3, 3);
  borated_pe_mat->AddElement(natured_b_el, 0.864);
  borated_pe_mat->AddElement(nist_manager->FindOrBuildElement("C"), 0.116);
  borated_pe_mat->AddElement(nist_manager->FindOrBuildElement("H"), 0.02);

  // Hydrogen Borated
  hb_mat = new G4Material("Hydrogen_Borated", 1.32 * g / cm3, 2);
  hb_mat->AddElement(natured_b_el, 0.804);
  hb_mat->AddElement(nist_manager->FindOrBuildElement("H"), 0.196); // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // define world
  G4Box* world_solid = new G4Box("World", 0.5 * WorldSizeX, 0.5 * WorldSizeY, 0.5 * WorldSizeZ);
  world_log = new G4LogicalVolume(world_solid, air_mat, "World");
  G4VPhysicalVolume* world_phys = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), world_log, "World", nullptr, false, 0, check_overlaps);

  // Reaction chamber geometry.
  // ChamberShell is the stainless-steel mechanical shell.
  // ChamberVacuum is the evacuated internal tracking volume, including
  // the central cavity and open beam/diagnostic ports. It is not a wall.
  G4LogicalVolume* chamber_shell_log = GetChamberShellLog("ChamberShell");
  G4LogicalVolume* chamber_vacuum_log = GetChamberVacuumLog("ChamberVacuum");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), chamber_shell_log, "ChamberShell", world_log, false, 0, check_overlaps);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), chamber_vacuum_log, "ChamberVacuum", world_log, false, 0, check_overlaps);

  // Target
  G4double z_target = TargetZPos;
  G4LogicalVolume* target_log = GetTargetLog("Target");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_target), target_log, "Target", chamber_vacuum_log, false, 0, check_overlaps);

  // Target Backing
  if (flag_target_backing) {
    G4double z_target_backing = z_target + target_thickness / 2. + TargetBackingThickness / 2.;
    G4LogicalVolume* target_backing_log = GetTargetBackingLog("TargetBacking");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, z_target_backing), target_backing_log, "TargetBacking", chamber_vacuum_log, false, 0, check_overlaps);

    G4Region* target_backing_reg = new G4Region("TargetBacking");
    target_backing_reg->AddRootLogicalVolume(target_backing_log);
  }

  // Non-physical scoring shell.  It is made from the same G4_Galactic
  // material as the chamber vacuum, so crossing it introduces no material
  // energy loss.  Its only role is to provide a well-defined sensitive
  // boundary outside the target and inside the silicon array.
  if (VirtualSphereConfig::GetEnabled()) {
    const G4double inner_radius = VirtualSphereConfig::GetRadius();
    const G4double outer_radius = inner_radius + VirtualSphereConfig::GetThickness();
    auto sphere_solid = new G4Sphere("VirtualSphereSolid", inner_radius, outer_radius,
                                     0. * deg, 360. * deg, 0. * deg, 180. * deg);
    virtual_sphere_log = new G4LogicalVolume(sphere_solid, vaccum_mat, "VirtualSphereLog");
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., TargetZPos), virtual_sphere_log,
                      "VirtualSphere", chamber_vacuum_log, false, 0, check_overlaps);

    auto sphere_vis = new G4VisAttributes(G4Colour(0.1, 0.8, 0.9, 0.25));
    sphere_vis->SetForceWireframe(true);
    virtual_sphere_log->SetVisAttributes(sphere_vis);
  }

  // G4Region for Cut
  G4Region* chamber_shell_reg = new G4Region("ChamberShell");
  chamber_shell_reg->AddRootLogicalVolume(chamber_shell_log);
  G4Region* target_reg = new G4Region("Target");
  target_reg->AddRootLogicalVolume(target_log);

  // Si Array
  si_array = new SiArray(chamber_vacuum_log);
  si_array->Construct();

  // HPGe Array
  hpge_array = new HPGeArray(world_log);
  hpge_array->Construct();

  // Labr3 Array
  labr3_array = new LaBr3Array(world_log);
  labr3_array->Construct();

  //
  G4RegionStore* reg_store = G4RegionStore::GetInstance();
  for (G4int i = 0; i < reg_store->size(); ++i) {
    G4Region* reg = (*reg_store)[i];
    G4cout << "G4Region for Cuts " << reg->GetName() << G4endl;
    if (!reg->GetProductionCuts()) {
      G4ProductionCuts* cuts = new G4ProductionCuts();
      cuts->SetProductionCut(1. * mm, "e-");
      cuts->SetProductionCut(1. * mm, "e+");
      reg->SetProductionCuts(cuts);
    }
  }

  //
  G4double max_step = StepMax4Proton;
  step_limit = new G4UserLimits(max_step);
  target_log->SetUserLimits(step_limit);

  return world_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  G4SDManager* sd_manager = G4SDManager::GetSDMpointer();

  auto si_sd = new SiSD("SiSD", "SiFrontHitCollection", "SiBackHitCollection");
  sd_manager->AddNewDetector(si_sd);
  if (si_array)
    si_array->MakeSensitive(si_sd);

  auto hpge_sd = new HPGeSD("HPGeSD", "HPGeHitCollection");
  sd_manager->AddNewDetector(hpge_sd);
  if (hpge_array)
    hpge_array->MakeSensitive(hpge_sd);

  auto labr3_sd = new LaBr3SD("LaBr3SD", "LaBr3HitCollection");
  sd_manager->AddNewDetector(labr3_sd);
  if (labr3_array)
    labr3_array->MakeSensitive(labr3_sd);

  if (virtual_sphere_log && VirtualSphereConfig::GetEnabled()) {
    auto virtual_sphere_sd = new VirtualSphereSD("VirtualSphereSD", "VirtualSphereHitCollection");
    sd_manager->AddNewDetector(virtual_sphere_sd);
    virtual_sphere_log->SetSensitiveDetector(virtual_sphere_sd);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetChamberShellLog(G4String name)
{
  G4cout << " --->: GetChamberShellLog, begin" << G4endl;

  // Build only the stainless-steel material volume of the chamber.
  // The evacuated inner space is constructed separately in
  // GetChamberVacuumLog() and placed as an adjacent daughter of World.

  G4double x_outer = ChamberX;
  G4double y_outer = ChamberY;
  G4double z_outer = ChamberZ;

  G4double x_inner = x_outer - 2. * ChamberThickness;
  G4double y_inner = y_outer - 2. * ChamberThickness;
  G4double z_inner = z_outer - 2. * ChamberThickness;

  G4double thickness = ChamberThickness;
  G4double flange_xr = FlangeXR;
  G4double flange_yr = FlangeYR;
  G4double flange_zr = FlangeZR;

  G4double flange_xxr = FlangeXR + thickness / 2.;
  G4double flange_yyr = FlangeYR + thickness / 2.;
  G4double flange_zzr = FlangeZR + thickness / 2.;

  G4double flange_xxr2 = FlangeXR2;
  G4double flange_yyr2 = FlangeYR2;
  G4double flange_zzr2 = FlangeZR2;

  G4double flange_xxrh = FlangeXRH;
  G4double flange_yyrh = FlangeXRH;
  G4double flange_zzrh = FlangeXRH;

  G4double rot_x_angle, rot_y_angle, rot_z_angle;

  //
  G4Box* outer_solid = new G4Box("outer_solid", x_outer / 2., y_outer / 2., z_outer / 2.);
  G4Box* inner_solid = new G4Box("inner_solid", x_inner / 2., y_inner / 2., z_inner / 2.);
  G4SubtractionSolid* step1_solid = new G4SubtractionSolid("step1_solid", outer_solid, inner_solid, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0))); // flange_x_solid
  G4cout << " --->: flange_x_solid, begin" << G4endl;
  G4Tubs* flange_x_solid = new G4Tubs("flange_x_solid", 0. * mm, flange_xr, thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_x_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 90. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_x_rot_matrix->rotateX(rot_x_angle);
  flange_x_rot_matrix->rotateY(rot_y_angle);
  flange_x_rot_matrix->rotateZ(rot_z_angle);
  flange_x_rot_matrix->print(G4cout);

  // 123.5 = (257-10)/2.
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_x_pos = G4ThreeVector(123.5 * mm, 0. * mm, 158.5 * mm);
  G4SubtractionSolid* step2_solid = new G4SubtractionSolid("step2_solid", step1_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(123.5 * mm, 0. * mm, -158.5 * mm);
  G4SubtractionSolid* step3_solid = new G4SubtractionSolid("step3_solid", step2_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(-123.5 * mm, 0. * mm, 158.5 * mm);
  G4SubtractionSolid* step4_solid = new G4SubtractionSolid("step4_solid", step3_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(-123.5 * mm, 0. * mm, -158.5 * mm);
  G4SubtractionSolid* step5_solid = new G4SubtractionSolid("step5_solid", step4_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos)); // flange_y_solid
  G4cout << " --->: flange_y_solid, begin" << G4endl;
  G4Tubs* flange_y_solid = new G4Tubs("flange_y_solid", 0. * mm, flange_yr, thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_y_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 90. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_y_rot_matrix->rotateX(rot_x_angle);
  flange_y_rot_matrix->rotateY(rot_y_angle);
  flange_y_rot_matrix->rotateZ(rot_z_angle);
  flange_y_rot_matrix->print(G4cout);

  // 123.5 = (257-10)/2.
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_y_pos = G4ThreeVector(0. * mm, 123.5 * mm, 158.5 * mm);
  G4SubtractionSolid* step6_solid = new G4SubtractionSolid("step6_solid", step5_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, 123.5 * mm, -158.5 * mm);
  G4SubtractionSolid* step7_solid = new G4SubtractionSolid("step7_solid", step6_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, -123.5 * mm, 158.5 * mm);
  G4SubtractionSolid* step8_solid = new G4SubtractionSolid("step8_solid", step7_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, -123.5 * mm, -158.5 * mm);
  G4SubtractionSolid* step9_solid = new G4SubtractionSolid("step9_solid", step8_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos)); // flange_z_solid
  G4cout << " --->: flange_z_solid, begin" << G4endl;
  G4Tubs* flange_z_solid = new G4Tubs("flange_z_solid", 0. * mm, flange_zr, thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_z_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_z_rot_matrix->rotateX(rot_x_angle);
  flange_z_rot_matrix->rotateY(rot_y_angle);
  flange_z_rot_matrix->rotateZ(rot_z_angle);
  flange_z_rot_matrix->print(G4cout);

  // 323.5 = (657-10)/2.
  G4ThreeVector flange_z_pos = G4ThreeVector(0. * mm, 0. * mm, 323.5 * mm);
  G4SubtractionSolid* step10_solid = new G4SubtractionSolid("step10_solid", step9_solid, flange_z_solid, G4Transform3D(*flange_z_rot_matrix, flange_z_pos));
  flange_z_pos = G4ThreeVector(0. * mm, 0. * mm, -323.5 * mm);
  G4SubtractionSolid* step11_solid = new G4SubtractionSolid("step11_solid", step10_solid, flange_z_solid, G4Transform3D(*flange_z_rot_matrix, flange_z_pos)); // flange_xx_solid
  G4cout << " --->: flange_xx_solid, begin" << G4endl;

  G4double xx_z_plane[6] = {-5., flange_xxrh, flange_xxrh, flange_xxrh + flange_xxrh - thickness, flange_xxrh + flange_xxrh - thickness, flange_xxrh + flange_xxrh};
  G4double xx_r_inner_plane[6] = {flange_xr, flange_xr, flange_xr, flange_xr, 0, 0};
  G4double xx_r_outer_plane[6] = {flange_xxr, flange_xxr, flange_xxr2, flange_xxr2, flange_xxr2, flange_xxr2};

  G4Polycone* flange_xx_solid = new G4Polycone("flange_xx_solid", 0. * deg, 360. * deg, 6, xx_z_plane, xx_r_inner_plane, xx_r_outer_plane);

  G4RotationMatrix* flange_xx_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 90. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_xx_rot_matrix->rotateX(rot_x_angle);
  flange_xx_rot_matrix->rotateY(rot_y_angle);
  flange_xx_rot_matrix->rotateZ(rot_z_angle);
  flange_xx_rot_matrix->print(G4cout);

  // 128.5 = 257/2.
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_xx_pos = G4ThreeVector(128.5 * mm, 0. * mm, 158.5 * mm);
  G4UnionSolid* step12_solid = new G4UnionSolid("step12_solid", step11_solid, flange_xx_solid, G4Transform3D(*flange_xx_rot_matrix, flange_xx_pos));
  flange_xx_pos = G4ThreeVector(128.5 * mm, 0. * mm, -158.5 * mm);
  G4UnionSolid* step13_solid = new G4UnionSolid("step13_solid", step12_solid, flange_xx_solid, G4Transform3D(*flange_xx_rot_matrix, flange_xx_pos));

  flange_xx_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = -90. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_xx_rot_matrix->rotateX(rot_x_angle);
  flange_xx_rot_matrix->rotateY(rot_y_angle);
  flange_xx_rot_matrix->rotateZ(rot_z_angle);
  flange_xx_rot_matrix->print(G4cout);
  flange_xx_pos = G4ThreeVector(-128.5 * mm, 0. * mm, -158.5 * mm);
  G4UnionSolid* step14_solid = new G4UnionSolid("step14_solid", step13_solid, flange_xx_solid, G4Transform3D(*flange_xx_rot_matrix, flange_xx_pos));
  flange_xx_pos = G4ThreeVector(-128.5 * mm, 0. * mm, 158.5 * mm);
  G4UnionSolid* step15_solid = new G4UnionSolid("step15_solid", step14_solid, flange_xx_solid, G4Transform3D(*flange_xx_rot_matrix, flange_xx_pos)); // flange_yy_solid
  G4cout << " --->: flange_yy_solid, begin" << G4endl;

  G4double yy_z_plane[6] = {-5., flange_yyrh, flange_yyrh, flange_yyrh + flange_yyrh - thickness, flange_yyrh + flange_yyrh - thickness, flange_yyrh + flange_yyrh};
  G4double yy_r_inner_plane[6] = {flange_yr, flange_yr, flange_yr, flange_yr, 0, 0};
  G4double yy_r_outer_plane[6] = {flange_yyr, flange_yyr, flange_yyr2, flange_yyr2, flange_yyr2, flange_yyr2};

  G4Polycone* flange_yy_solid = new G4Polycone("flange_yy_solid", 0. * deg, 360. * deg, 6, yy_z_plane, yy_r_inner_plane, yy_r_outer_plane);

  G4RotationMatrix* flange_yy_rot_matrix = new G4RotationMatrix();
  rot_x_angle = -90. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_yy_rot_matrix->rotateX(rot_x_angle);
  flange_yy_rot_matrix->rotateY(rot_y_angle);
  flange_yy_rot_matrix->rotateZ(rot_z_angle);
  flange_yy_rot_matrix->print(G4cout);

  // 128.5 = 257/2.
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_yy_pos = G4ThreeVector(0. * mm, 128.5 * mm, 158.5 * mm);
  G4UnionSolid* step16_solid = new G4UnionSolid("step16_solid", step15_solid, flange_yy_solid, G4Transform3D(*flange_yy_rot_matrix, flange_yy_pos));
  flange_yy_pos = G4ThreeVector(0. * mm, 128.5 * mm, -158.5 * mm);
  G4UnionSolid* step17_solid = new G4UnionSolid("step17_solid", step16_solid, flange_yy_solid, G4Transform3D(*flange_yy_rot_matrix, flange_yy_pos));

  flange_yy_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 90. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_yy_rot_matrix->rotateX(rot_x_angle);
  flange_yy_rot_matrix->rotateY(rot_y_angle);
  flange_yy_rot_matrix->rotateZ(rot_z_angle);
  flange_yy_rot_matrix->print(G4cout);
  flange_yy_pos = G4ThreeVector(0. * mm, -128.5 * mm, -158.5 * mm);
  G4UnionSolid* step18_solid = new G4UnionSolid("step18_solid", step17_solid, flange_yy_solid, G4Transform3D(*flange_yy_rot_matrix, flange_yy_pos));
  flange_yy_pos = G4ThreeVector(0. * mm, -128.5 * mm, 158.5 * mm);
  G4UnionSolid* step19_solid = new G4UnionSolid("step19_solid", step18_solid, flange_yy_solid, G4Transform3D(*flange_yy_rot_matrix, flange_yy_pos)); // flange_zz_solid
  G4cout << " --->: flange_zz_solid, begin" << G4endl;

  G4double zz_z_plane[6] = {-5., flange_zzrh, flange_zzrh, flange_zzrh + flange_zzrh - thickness, flange_zzrh + flange_zzrh - thickness, flange_zzrh + flange_zzrh};
  G4double zz_r_inner_plane[6] = {flange_zr, flange_zr, flange_zr, flange_zr, 0, 0};
  G4double zz_r_outer_plane[6] = {flange_zzr, flange_zzr, flange_zzr2, flange_zzr2, flange_zzr2, flange_zzr2};
  G4cout << " --->: flange_zz_solid, begin" << G4endl;

  G4Polycone* flange_zz_solid = new G4Polycone("flange_zz_solid", 0. * deg, 360. * deg, 6, zz_z_plane, zz_r_inner_plane, zz_r_outer_plane);

  G4RotationMatrix* flange_zz_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_zz_rot_matrix->rotateX(rot_x_angle);
  flange_zz_rot_matrix->rotateY(rot_y_angle);
  flange_zz_rot_matrix->rotateZ(rot_z_angle);
  flange_zz_rot_matrix->print(G4cout);

  // 328.5 = 657/2.
  G4ThreeVector flange_zz_pos = G4ThreeVector(0. * mm, 0. * mm, 328.5 * mm);
  G4UnionSolid* step20_solid = new G4UnionSolid("step20_solid", step19_solid, flange_zz_solid, G4Transform3D(*flange_zz_rot_matrix, flange_zz_pos));

  flange_zz_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 180. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_zz_rot_matrix->rotateX(rot_x_angle);
  flange_zz_rot_matrix->rotateY(rot_y_angle);
  flange_zz_rot_matrix->rotateZ(rot_z_angle);
  flange_zz_rot_matrix->print(G4cout);
  flange_zz_pos = G4ThreeVector(0. * mm, 0. * mm, -328.5 * mm);
  G4UnionSolid* step21_solid = new G4UnionSolid("step21_solid", step20_solid, flange_zz_solid, G4Transform3D(*flange_zz_rot_matrix, flange_zz_pos));

  TString log_name = TString::Format("%s_log", name.c_str());
  G4LogicalVolume* log = new G4LogicalVolume(step21_solid, stainless_steel_mat, log_name.Data());

  // color
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.3, 0.4, 0.5, 0.8));
  vis_att->SetForceSolid(true);
  log->SetVisAttributes(vis_att);

  return log;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetChamberVacuumLog(G4String name)
{
  G4cout << " --->: GetChamberVacuumLog, begin" << G4endl;

  // Build the evacuated internal volume of the reaction chamber.
  // This logical volume represents the beam/particle tracking space,
  // not an additional wall layer. Target and internal Si detectors are
  // placed inside this vacuum volume.

  G4double x = ChamberX - 2. * ChamberThickness;
  G4double y = ChamberY - 2. * ChamberThickness;
  G4double z = ChamberZ - 2. * ChamberThickness;

  G4double x_thickness = ChamberThickness + 2. * FlangeXRH;
  G4double y_thickness = ChamberThickness + 2. * FlangeYRH;
  G4double z_thickness = ChamberThickness + 2. * FlangeZRH;
  G4double flange_xr = FlangeXR;
  G4double flange_yr = FlangeYR;
  G4double flange_zr = FlangeZR;

  G4double rot_x_angle, rot_y_angle, rot_z_angle;

  //
  G4Box* step1_solid = new G4Box("step1_solid", x / 2., y / 2., z / 2.); // flange_x_solid
  G4cout << " --->: flange_x_solid, begin" << G4endl;
  G4Tubs* flange_x_solid = new G4Tubs("flange_x_solid", 0. * mm, flange_xr, x_thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_x_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 90. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_x_rot_matrix->rotateX(rot_x_angle);
  flange_x_rot_matrix->rotateY(rot_y_angle);
  flange_x_rot_matrix->rotateZ(rot_z_angle);
  flange_x_rot_matrix->print(G4cout);

  // 103.5 = (257-10)/2. - 10 - 10
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_x_pos = G4ThreeVector(103.5 * mm, 0. * mm, 158.5 * mm);
  G4UnionSolid* step2_solid = new G4UnionSolid("step2_solid", step1_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(103.5 * mm, 0. * mm, -158.5 * mm);
  G4UnionSolid* step3_solid = new G4UnionSolid("step3_solid", step2_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(-103.5 * mm, 0. * mm, 158.5 * mm);
  G4UnionSolid* step4_solid = new G4UnionSolid("step4_solid", step3_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos));
  flange_x_pos = G4ThreeVector(-103.5 * mm, 0. * mm, -158.5 * mm);
  G4UnionSolid* step5_solid = new G4UnionSolid("step5_solid", step4_solid, flange_x_solid, G4Transform3D(*flange_x_rot_matrix, flange_x_pos)); // flange_y_solid
  G4cout << " --->: flange_y_solid, begin" << G4endl;
  G4Tubs* flange_y_solid = new G4Tubs("flange_y_solid", 0. * mm, flange_yr, y_thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_y_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 90. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_y_rot_matrix->rotateX(rot_x_angle);
  flange_y_rot_matrix->rotateY(rot_y_angle);
  flange_y_rot_matrix->rotateZ(rot_z_angle);
  flange_y_rot_matrix->print(G4cout);

  // 103.5 = (257-10)/2. -10. - 10.
  // 158.5 = 657/2. - 170
  G4ThreeVector flange_y_pos = G4ThreeVector(0. * mm, 103.5 * mm, 158.5 * mm);
  G4UnionSolid* step6_solid = new G4UnionSolid("step6_solid", step5_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, 103.5 * mm, -158.5 * mm);
  G4UnionSolid* step7_solid = new G4UnionSolid("step7_solid", step6_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, -103.5 * mm, 158.5 * mm);
  G4UnionSolid* step8_solid = new G4UnionSolid("step8_solid", step7_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos));
  flange_y_pos = G4ThreeVector(0. * mm, -103.5 * mm, -158.5 * mm);
  G4UnionSolid* step9_solid = new G4UnionSolid("step9_solid", step8_solid, flange_y_solid, G4Transform3D(*flange_y_rot_matrix, flange_y_pos)); // flange_z_solid
  G4cout << " --->: flange_z_solid, begin" << G4endl;
  G4Tubs* flange_z_solid = new G4Tubs("flange_z_solid", 0. * mm, flange_zr, z_thickness, 0. * deg, 360. * deg);

  G4RotationMatrix* flange_z_rot_matrix = new G4RotationMatrix();
  rot_x_angle = 0. * deg;
  rot_y_angle = 0. * deg;
  rot_z_angle = 0. * deg;
  G4cout << "x angle " << rot_x_angle << "y angle " << rot_y_angle << " z angle " << rot_z_angle << G4endl;
  flange_z_rot_matrix->rotateX(rot_x_angle);
  flange_z_rot_matrix->rotateY(rot_y_angle);
  flange_z_rot_matrix->rotateZ(rot_z_angle);
  flange_z_rot_matrix->print(G4cout);

  // 303.5 = (657-10)/2. - 10. - 10.
  G4ThreeVector flange_z_pos = G4ThreeVector(0. * mm, 0. * mm, 303.5 * mm);
  G4UnionSolid* step10_solid = new G4UnionSolid("step10_solid", step9_solid, flange_z_solid, G4Transform3D(*flange_z_rot_matrix, flange_z_pos));
  flange_z_pos = G4ThreeVector(0. * mm, 0. * mm, -303.5 * mm);
  G4UnionSolid* step11_solid = new G4UnionSolid("step11_solid", step10_solid, flange_z_solid, G4Transform3D(*flange_z_rot_matrix, flange_z_pos));

  TString log_name = TString::Format("%s_log", name.c_str());
  G4LogicalVolume* log = new G4LogicalVolume(step11_solid, vaccum_mat, log_name.Data());

  // color
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.4, 0.3, 0.6));
  vis_att->SetForceSolid(true);
  log->SetVisAttributes(vis_att);

  return log;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetTargetLog(G4String name)
{
  G4cout << " --->: GetTargetLog, begin" << G4endl;

  G4double r = TargetR;

  G4Tubs* solid = new G4Tubs("solid", 0. * mm, r, target_thickness / 2., 0. * deg, 360. * deg);

  TString log_name = TString::Format("%s_log", name.c_str());
  G4LogicalVolume* log = new G4LogicalVolume(solid, target_mat, log_name.Data());

  // color
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.5, 0.6, 0.3, 0.8));
  vis_att->SetForceSolid(true);
  log->SetVisAttributes(vis_att);

  return log;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* DetectorConstruction::GetTargetBackingLog(G4String name)
{
  G4cout << " --->: GetTargetBackingLog, begin" << G4endl;

  G4double r = TargetBackingR;
  G4double thickness = TargetBackingThickness;

  G4Tubs* solid = new G4Tubs("solid", 0. * mm, r, thickness / 2., 0. * deg, 360. * deg);

  TString log_name = TString::Format("%s_log", name.c_str());
  G4LogicalVolume* log = new G4LogicalVolume(solid, si_mat, log_name.Data());

  // color
  G4VisAttributes* vis_att = new G4VisAttributes(G4Colour(0.8, 0.8, 0.3, 0.6));
  vis_att->SetForceSolid(true);
  log->SetVisAttributes(vis_att);

  return log;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMaxStep(G4double max_step)
{
  if ((step_limit) && (max_step > 0.))
    step_limit->SetMaxAllowedStep(max_step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCheckOverlaps(G4bool co)
{
  check_overlaps = co;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTargetMaterial(G4String str)
{
  G4Material* mat = G4Material::GetMaterial(str);
  if (!mat) {
    G4cout << "cannot find " << str << " material " << G4endl;
    return;
  }
  target_mat = mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// /target/ messenger command wrappers.  Thickness / material / backing all
// change the built geometry, so they are only meaningful before the geometry
// is constructed (/run/initialize) or when followed by
// /run/reinitializeGeometry.  We flag the geometry as modified so that a
// subsequent /run/beamOn triggers a rebuild automatically.
void DetectorConstruction::SetTargetThicknessCmd(G4double th)
{
  if (th <= 0.) {
    G4cout << "[Target] ignoring non-positive thickness " << th / um << " um" << G4endl;
    return;
  }
  SetTargetThickness(th);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "[Target] thickness set to " << target_thickness / um << " um"
         << " (geometry flagged for rebuild)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTargetMaterialCmd(G4String str)
{
  SetTargetMaterial(str);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "[Target] material set to " << str
         << " (geometry flagged for rebuild)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTargetBackingFlagCmd(G4bool bl)
{
  SetTargetBackingFlag(bl);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "[Target] backing " << (bl ? "enabled" : "disabled")
         << " (geometry flagged for rebuild)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintTargetConfigCommand()
{
  G4cout << "==== Target configuration ====" << G4endl
         << "  thickness = " << target_thickness / um << " um"
         << " (" << target_thickness / mm << " mm)" << G4endl
         << "  material  = "
         << (target_mat ? target_mat->GetName() : G4String("<null>")) << G4endl
         << "  backing   = " << (flag_target_backing ? "on" : "off") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineCommands()
{
  target_messenger = std::make_unique<G4GenericMessenger>(
      this, "/target/", "Target geometry controls (thickness / material / backing)");

  target_messenger->DeclareMethodWithUnit("thickness", "um",
                                          &DetectorConstruction::SetTargetThicknessCmd,
                                          "Target thickness (geometry; set before /run/initialize "
                                          "or follow with /run/reinitializeGeometry). "
                                          "Scan this for triple-coincidence efficiency vs thickness.");
  target_messenger->DeclareMethod("material",
                                  &DetectorConstruction::SetTargetMaterialCmd,
                                  "Target material name (geometry; must already be defined)");
  target_messenger->DeclareMethod("backing",
                                  &DetectorConstruction::SetTargetBackingFlagCmd,
                                  "Enable/disable the target backing (geometry)");
  target_messenger->DeclareMethod("printConfig",
                                  &DetectorConstruction::PrintTargetConfigCommand,
                                  "Print current target configuration");
}
