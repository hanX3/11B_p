#ifndef SiDetector_H
#define SiDetector_H 1

#include "Constants.hh"
#include "DetectorChannel.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <array>
#include <map>

class SiDetector
{
public:
  SiDetector(G4LogicalVolume* log);
  ~SiDetector();

  void SetName(const G4String name) { si_name = name; }
  void SetRingId(const G4int id) { ring_id = id; }
  void SetSectorId(const G4int id) { sector_id = id; }

  void ConstructSiDetector(const std::array<G4double, 2>& si_par, G4Material* mat);
  void ConstructSiBoxDetector(const std::array<G4double, 3>& si_box_par, G4Material* mat);
  void ConstructSiWedgeDetector(const std::array<G4double, 4>& si_wedge_par, G4Material* mat);
  void PlaceSiDetector(G4RotationMatrix* rot, const G4ThreeVector& pos);
  void PlaceSiDetector(const G4Transform3D& transform_3d);


  G4String GetName() { return si_name; }
  G4int GetRingId() { return ring_id; }
  G4int GetSectorId() { return sector_id; }
  G4LogicalVolume* GetLog() { return si_detector_log; }
  G4VPhysicalVolume* GetPhy() { return si_detector_phy; }

  G4LogicalVolume* exp_hall_log;

private:
  G4LogicalVolume* si_detector_log = nullptr;
  G4VPhysicalVolume* si_detector_phy = nullptr;

  G4String si_name;
  G4int ring_id = -1;
  G4int sector_id = -1;
  G4bool check_overlaps = true;

public:
  static std::map<G4String, G4int> map_name_to_ring_id;
  static std::map<G4int, G4String> map_ring_id_to_name;
  static std::map<G4String, G4int> map_name_to_sectors;
  static std::map<G4String, std::array<G4double, 2>> map_si_par;
  static std::map<G4String, std::array<G4double, 3>> map_si_box_par;
  // inner full width, outer full width, active length, thickness [mm]
  static std::map<G4String, std::array<G4double, 4>> map_si_wedge_par;
  static std::map<G4String, G4double> map_si_inner_radius;
  static std::map<G4String, std::array<G4double, 3>> map_placement_par;
  static std::map<G4String, std::array<G4double, 4>> map_color_par;
};

#endif
