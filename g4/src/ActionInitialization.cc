#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootIO.hh"
#include "RunAction.hh"
#include "TrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ActionInitialization::ActionInitialization(ULong64_t seed)
    : random_seed(seed) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ActionInitialization::~ActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(nullptr));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ActionInitialization::Build() const
{
  auto root_io = new RootIO();
  root_io->SetRandomSeed(random_seed);

  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction(root_io));
  SetUserAction(new EventAction(root_io));
  SetUserAction(new TrackingAction());
}
