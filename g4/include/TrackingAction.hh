#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "Constants.hh"
#include "DataStructure.hh"
#include "G4UserTrackingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;
class RootIO;

//
class TrackingAction : public G4UserTrackingAction {
public:
  TrackingAction(DetectorConstruction* dc, RootIO* rio);
  ~TrackingAction() {};

  void  PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);

private:
  DetectorConstruction* detector;
  RootIO* root_io;
  TrackData track_data;
};

#endif

