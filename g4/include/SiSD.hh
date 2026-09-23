#ifndef SiSD_h
#define SiSD_h 1

#include "Constants.hh"
#include "SiHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"

#include <vector>

//
class SiSD : public G4VSensitiveDetector
{
public:
  SiSD(const G4String &name, const G4String &hits_collection_name);
  ~SiSD() override = default;

  // methods from base class
  void Initialize(G4HCofThisEvent *) override;
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
  void EndOfEvent(G4HCofThisEvent *) override;

private:
  SiHitsCollection *hits_collection;
  G4int hc_id;
};


#endif
