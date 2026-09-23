#ifndef LaBr3Array_h
#define LaBr3Array_h 1

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "Constants.hh"
#include "LaBr3Detector.hh"
#include "LaBr3SD.hh"

#include <vector>

class LaBr3Array
{
public:
  LaBr3Array(G4LogicalVolume *log);
  ~LaBr3Array();

public:
  void Construct();
  void MakeSensitive(LaBr3SD *labr3_sd);

public:
  void PrintDetectorDimensionInfo();

public:
  G4LogicalVolume *exp_hall_log;

private:
  G4Transform3D CalculatePlacement(G4String name, G4int sector_id);

private:
  G4Material *labr3_mat;
  G4Material *al_mat;

private:
  G4int labr3_numbers;

  std::vector<LaBr3Detector*> v_labr3_detector;
};

#endif