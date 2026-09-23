#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "Constants.hh"
#include "LaBr3SD.hh"
#include "LaBr3Array.hh"
#include "HPGeSD.hh"
#include "HPGeArray.hh"
#include "SiSD.hh"
#include "SiArray.hh"

#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

#include "G4Element.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"

#include <memory>

class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

public:
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

public:
  G4LogicalVolume* GetChamberShellLog(G4String name);
  G4LogicalVolume* GetChamberVacuumLog(G4String name);
  G4LogicalVolume* GetTargetLog(G4String name);
  G4LogicalVolume* GetTargetBackingLog(G4String name);

  G4LogicalVolume* GetFlangeYLog(G4String name);
  G4LogicalVolume* GetFlangeZLog(G4String name);

  // Set methods
  void SetMaxStep(G4double);
  void SetCheckOverlaps(G4bool);

  void SetTargetThickness(G4double th)
  {
    target_thickness = th;
  }
  G4double GetTargetThickness() const { return target_thickness; }
  void SetTargetMaterial(G4String str);

  void SetTargetBackingFlag(G4bool bl)
  {
    flag_target_backing = bl;
  }

  // ---- /target/ messenger command wrappers ----
  // thickness is geometry: set it BEFORE /run/initialize, or set it and issue
  // /run/reinitializeGeometry.  material and backing likewise affect the build.
  void SetTargetThicknessCmd(G4double th);
  void SetTargetMaterialCmd(G4String str);
  void SetTargetBackingFlagCmd(G4bool bl);
  void PrintTargetConfigCommand();

private:
  // methods
  void DefineMaterials();
  void DefineCommands();
  G4VPhysicalVolume* DefineVolumes();

  std::unique_ptr<G4GenericMessenger> target_messenger;

private:
  SiArray* si_array;
  HPGeArray* hpge_array;
  LaBr3Array* labr3_array;

private:
  G4bool flag_target_backing;

  G4double target_thickness;
  G4Material* target_mat;

private:
  //
  G4LogicalVolume* world_log;
  G4LogicalVolume* virtual_sphere_log;

  //
  G4Material* air_mat;
  G4Material* vaccum_mat;
  G4Material* al_mat;
  G4Material* si_mat;
  G4Material* stainless_steel_mat;
  G4Material* enriched_11b_mat;
  G4Material* natured_11b_high_density_mat;
  G4Material* natured_11b_low_density_mat;
  G4Material* borated_pe_mat;
  G4Material* hb_mat;

  G4UserLimits* step_limit;
  G4bool check_overlaps;
};

#endif
