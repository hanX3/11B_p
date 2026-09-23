#ifndef B2EventAction_H
#define B2EventAction_H 1

#include "Constants.hh"
#include "DataStructure.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class RootIO;
class TFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RootIO* root_io);
  ~EventAction() override;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

private:
  void GausEnergy(G4double res);
  bool IfThresholdTrigger(G4double threshold);

private:
  G4int hc_id_si;
  G4int hc_id_hpge;
  G4int hc_id_labr3;

  G4double threshold_si;
  G4double energy_resolution_si;
  G4double threshold_hpge;
  G4double energy_resolution_hpge;
  G4double threshold_labr3;
  G4double energy_resolution_labr3;

private:
  EventData event_data;

  RootIO* root_io;
};

#endif
