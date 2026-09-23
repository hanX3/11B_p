#ifndef HPGeSD_h
#define HPGeSD_h 1

#include "Constants.hh"
#include "HPGeHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"

#include <vector>

//
class HPGeSD : public G4VSensitiveDetector
{
public:
  HPGeSD(const G4String &name, const G4String &hits_collection_name);
  ~HPGeSD() override = default;

  // methods from base class
  void Initialize(G4HCofThisEvent *) override;
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
  void EndOfEvent(G4HCofThisEvent *) override;

private:
  HPGeHitsCollection *hits_collection;
  G4int hc_id;
};


#endif
