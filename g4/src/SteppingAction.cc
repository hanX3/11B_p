#include "SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step*)
{
  // Intentionally empty. Target proton step limiting is handled by
  // G4UserLimits together with ProtonStepLimiterPhysics.
}
