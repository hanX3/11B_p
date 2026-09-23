#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "DataStructure.hh"

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4LogicalVolume;

class RootIO;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(RootIO* root_io);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step*);

private:
  RootIO* root_io;
  StepData step_data;
};

#endif
