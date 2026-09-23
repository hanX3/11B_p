#ifndef SiArray_H
#define SiArray_H 1

#include "Constants.hh"
#include "DataStructure.hh"
#include "SiDetector.hh"
#include "SiSD.hh"

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"

#include <vector>

class SiArray
{
public:
  SiArray(G4LogicalVolume* log);
  ~SiArray();

public:
  void Construct();
  void MakeSensitive(SiSD* si_sd);

  // Build the static mapping from a front/back strip pair to the centre of the
  // corresponding ideal DSSD pixel.  The map contains geometry only and is
  // regenerated at run start using the current strip-count configuration.
  static std::vector<SiPixelMapData> BuildPixelMap();

public:
  void PrintDetectorDimensionInfo();

public:
  G4LogicalVolume* exp_hall_log;

private:
  struct DetectorGeometry
  {
    G4int detector_id = -1;
    G4int detector_model = 0; // 1: W1, 2: S3
    G4int subarray_id = -1;
    G4int module_id = -1;
    G4Transform3D local_to_world;
    G4double half_x = 0.;
    G4double half_y = 0.;
    G4double r_inner = 0.;
    G4double r_outer = 0.;
  };

  G4Transform3D CalculatePlacement(G4String name, G4int sector_id);
  void RegisterDetectorGeometry(SiDetector& detector, const G4Transform3D& transform);

private:
  G4Material* si_mat;
  G4Material* al_mat;

private:
  G4int si_numbers;

  std::vector<SiDetector*> v_si_detector;
  static std::vector<DetectorGeometry> detector_geometry;
};

#endif
