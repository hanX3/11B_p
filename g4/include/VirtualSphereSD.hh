#ifndef VirtualSphereSD_H
#define VirtualSphereSD_H 1

#include "VirtualSphereHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VSensitiveDetector.hh"

#include <unordered_set>

class VirtualSphereSD : public G4VSensitiveDetector
{
public:
  VirtualSphereSD(const G4String& name, const G4String& collection_name);
  ~VirtualSphereSD() override = default;

  void Initialize(G4HCofThisEvent*) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;
  void EndOfEvent(G4HCofThisEvent*) override;

private:
  VirtualSphereHitsCollection* hits_collection = nullptr;
  G4int hc_id = -1;
  std::unordered_set<G4int> recorded_track_ids;
};

#endif
