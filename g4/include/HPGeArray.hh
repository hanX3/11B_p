#ifndef HPGeArray_h
#define HPGeArray_h 1

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "Constants.hh"
#include "HPGeDetector.hh"
#include "HPGeSD.hh"

#include <vector>

class HPGeArray
{
public:
  HPGeArray(G4LogicalVolume *log);
  ~HPGeArray();

public:
  void Construct();
  void MakeSensitive(HPGeSD *hpge_sd);

public:
  void PrintDetectorDimensionInfo();

public:
  G4LogicalVolume *exp_hall_log;

private:
  G4Transform3D CalculatePlacement(G4String name, G4int sector_id);

private:
  G4Material *hpge_mat;
  G4Material *al_mat;

private:
  G4int hpge_numbers;

  std::vector<HPGeDetector*> v_hpge_detector;
};

#endif
