#ifndef LaBr3Detector_h
#define LaBr3Detector_h 1

#include "Constants.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4SystemOfUnits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class LaBr3Detector
{
public:
  LaBr3Detector(G4LogicalVolume *log);
  ~LaBr3Detector();

public:
  void SetName(const G4String name) { labr3_name = name; }
  void SetRingId(const G4int id) { ring_id = id; }
  void SetSectorId(const G4int id) { sector_id = id; }
  
  // LaBr3
  void ConstructLaBr3Detector(const std::array<G4double, 2> &labr3_par, G4Material *mat);
  void PlaceLaBr3Detector(G4RotationMatrix *rot, const G4ThreeVector &pos);
  void PlaceLaBr3Detector(const G4Transform3D &transfrom_3d);

  // Al Shell
  void ConstructAlShell(const std::array<G4double, 3> &al_par, G4Material *mat);
  void PlaceAlShell(G4RotationMatrix *rot, const G4ThreeVector &pos);
  void PlaceAlShell(const G4Transform3D &transfrom_3d);

public:
  G4String GetName() { return labr3_name; }
  G4int GetRingId() { return ring_id; }
  G4int GetSectorId() { return sector_id; }

  G4LogicalVolume *GetLog() { return labr3_detector_log; }
  G4VPhysicalVolume *GetPhy() { return labr3_detector_phy; }

public:
  G4LogicalVolume *exp_hall_log;

private:
  G4LogicalVolume *labr3_detector_log;
  G4VPhysicalVolume *labr3_detector_phy;

  G4LogicalVolume *al_shell_log;
  G4VPhysicalVolume *al_shell_phy;

private:
  G4String labr3_name;
  G4int ring_id; // 1,
  G4int sector_id; // 1,

private:
  G4bool check_overlaps;

public:
  static std::map<G4String, G4int> map_name_to_ring_id;
  static std::map<G4int, G4String> map_ring_id_to_name;
  static std::map<G4String, G4int> map_name_to_sectors;
  static std::map<G4String, std::array<G4double, 2>> map_labr3_par;
  static std::map<G4String, std::array<G4double, 3>> map_placement_par;
  static std::map<G4String, std::array<G4double, 3>> map_al_par;
  static std::map<G4String, std::array<G4double, 4>> map_color_par;
  static std::map<G4String, std::array<G4double, 4>> map_color_al_shell_par;
};


#endif