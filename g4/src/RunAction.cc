#include "RunAction.hh"
#include "RootIO.hh"
#include "OutputConfig.hh"
#include "VirtualSphereConfig.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4Timer.hh"
#include "unistd.h"
#include <fstream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction(RootIO* r_io)
    : G4UserRunAction(), root_io(r_io)
{
  // set printing event number per each 100000 events
  G4RunManager::GetRunManager()->SetPrintProgress(100000);

  timer = new G4Timer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::~RunAction()
{
  delete timer;
  timer = nullptr;

  delete root_io;
  root_io = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::BeginOfRunAction(const G4Run* run)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  if (root_io) {
    if (OutputConfig::GetSaveReaction() || (OutputConfig::GetSaveVirtualSphere() && VirtualSphereConfig::GetEnabled()) || OutputConfig::GetSaveEvent()) {
      root_io->OpenDataFile();
      G4cout << "open unified event ROOT file" << G4endl;
    }
  }

  int run_id = run->GetRunID();
  timer->Start();
  G4cout << "======================   RunID = " << run_id << "  ======================" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run* run)
{
  if (root_io) {
    if (OutputConfig::GetSaveReaction() || (OutputConfig::GetSaveVirtualSphere() && VirtualSphereConfig::GetEnabled()) || OutputConfig::GetSaveEvent()) {
      root_io->CloseDataFile();
    }
  }

  // Print results
  G4cout << "  The run was " << run->GetNumberOfEvent() << " events " << G4endl;

  timer->Stop();
  G4cout << " time:  " << *timer << G4endl;
}
