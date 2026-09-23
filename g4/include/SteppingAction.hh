#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "G4UserSteppingAction.hh"

class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Reserved project hook for future per-step logic.
// The proton step limit in the target is implemented by G4UserLimits and
// ProtonStepLimiterPhysics; it does not depend on this user action.
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction() = default;
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step*) override;
};

#endif
