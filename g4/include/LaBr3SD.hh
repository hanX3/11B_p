#ifndef LaBr3SD_h
#define LaBr3SD_h 1

#include "Constants.hh"
#include "LaBr3Hit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"

#include <vector>

//
class LaBr3SD : public G4VSensitiveDetector
{
public:
  LaBr3SD(const G4String &name, const G4String &hits_collection_name);
  ~LaBr3SD() override = default;

  // methods from base class
  void Initialize(G4HCofThisEvent *) override;
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
  void EndOfEvent(G4HCofThisEvent *) override;

private:
  LaBr3HitsCollection *hits_collection;
  G4int hc_id;
};


#endif
