#ifndef SiSD_H
#define SiSD_H 1

#include "SiHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VSensitiveDetector.hh"

#include <map>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SiSD : public G4VSensitiveDetector
{
public:
  SiSD(const G4String& name, const G4String& front_collection_name, const G4String& back_collection_name);
  ~SiSD() override = default;

  void Initialize(G4HCofThisEvent*) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
  void EndOfEvent(G4HCofThisEvent*) override;

private:
  using StripKey = std::pair<G4int, G4int>;   // detector_id, strip_id
  using DepositKey = std::pair<G4int, G4int>; // detector_id, track_id

  void AccumulateStripHit(SiHitsCollection* collection,
                          std::map<StripKey, G4int>& hit_index,
                          const StripKey& key,
                          G4int detector_id,
                          SiReadoutSide side,
                          G4int strip_id,
                          G4double e_dep,
                          G4Step* step);

  void AccumulateDepositHit(const DepositKey& key,
                            G4int detector_id,
                            G4double e_dep,
                            const G4ThreeVector& entry_position,
                            const G4ThreeVector& deposit_position,
                            G4Step* step);

private:
  SiHitsCollection* front_hits_collection = nullptr;
  SiHitsCollection* back_hits_collection = nullptr;
  SiHitsCollection* deposit_hits_collection = nullptr;
  G4int front_hc_id = -1;
  G4int back_hc_id = -1;
  G4int deposit_hc_id = -1;

  // A physical energy deposit induces one signal on each DSSD face.  Within an
  // event, deposits sharing the same detector and strip are summed independently
  // in the front and back collections, preserving real strip-level ambiguity.
  std::map<StripKey, G4int> front_hit_index;
  std::map<StripKey, G4int> back_hit_index;

  // One physical entry per Geant4 track in one physical Si detector.
  std::map<DepositKey, G4int> deposit_hit_index;
};

#endif
