#ifndef TrackingAction_H
#define TrackingAction_H 1

#include "Constants.hh"
#include "DataStructure.hh"
#include "G4UserTrackingAction.hh"
#include "globals.hh"

class RootIO;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class TrackingAction : public G4UserTrackingAction
{
public:
  explicit TrackingAction(RootIO* root_io);
  ~TrackingAction(){};

  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);

private:
  RootIO* root_io;
  TrackData track_data;
};

#endif
