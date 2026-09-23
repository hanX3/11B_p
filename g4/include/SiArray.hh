#ifndef SiArray_h
#define SiArray_h 1

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "Constants.hh"
#include "SiDetector.hh"
#include "SiSD.hh"

#include <vector>

class SiArray
{
public:
  SiArray(G4LogicalVolume *log);
  ~SiArray();

public:
  void Construct();
  void MakeSensitive(SiSD *si_sd);

public:
  void PrintDetectorDimensionInfo();

public:
  G4LogicalVolume *exp_hall_log;

private:
  G4Transform3D CalculatePlacement(G4String name, G4int sector_id);

private:
  G4Material *si_mat;
  G4Material *al_mat;

private:
  G4int si_numbers;

  std::vector<SiDetector*> v_si_detector;
};

#endif
