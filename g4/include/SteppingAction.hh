#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "DataStructure.hh"

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;


class DetectorConstruction;
class EventAction;
class RootIO;

//
class SteppingAction: public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction* dc, RootIO* rio);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step*);

private:
  const DetectorConstruction* detector_construction;
  RootIO* root_io;
  StepData step_data;
};


#endif
