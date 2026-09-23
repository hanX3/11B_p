#ifndef B2EventAction_H
#define B2EventAction_H 1

#include "DataStructure.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

class RootIO;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RootIO* root_io);
  ~EventAction() override;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

private:
  G4int hc_id_si_front;
  G4int hc_id_si_back;
  G4int hc_id_si_deposit;
  G4int hc_id_hpge;
  G4int hc_id_labr3;
  G4int hc_id_virtual_sphere;

  EventData event_data;
  RootIO* root_io;
};

#endif
