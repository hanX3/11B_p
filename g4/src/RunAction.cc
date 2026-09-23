#include "RunAction.hh"
#include "RootIO.hh"
#include "PrimaryGeneratorAction.hh"
#include "OutputConfig.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleGun.hh"
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
    if (OutputConfig::GetSaveReaction()) {
      root_io->OpenReactionFile();
      G4cout << "open reaction root file " << G4endl;
    }
    if (OutputConfig::GetSaveEvent()) {
      root_io->OpenEventFile();
      G4cout << "open event root file " << G4endl;
    }
    if (OutputConfig::GetSaveTrack()) {
      root_io->OpenTrackFile();
      G4cout << "open track root file " << G4endl;
    }
    if (OutputConfig::GetSaveStep()) {
      root_io->OpenStepFile();
      G4cout << "open step root file " << G4endl;
    }
  }

  int run_id = run->GetRunID();
  timer->Start();
  G4cout << "======================   RunID = " << run_id << "  ======================" << G4endl;
  G4cout << "1 energy " << G4endl;
  G4cout << "2 track id " << G4endl;
  G4cout << "3 particle name " << G4endl;
  G4cout << "4 pre postion " << G4endl;
  G4cout << "5 post postion " << G4endl;
  G4cout << "6 detector name " << G4endl;
  G4cout << "7 step length " << G4endl;
  G4cout << "8 delat postion " << G4endl;
  G4cout << "9 delat time " << G4endl;
  G4cout << "10 delat energy " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::EndOfRunAction(const G4Run* run)
{
  if (root_io) {
    if (OutputConfig::GetSaveReaction()) {
      root_io->CloseReactionFile();
    }
    if (OutputConfig::GetSaveEvent()) {
      root_io->CloseEventFile();
    }
    if (OutputConfig::GetSaveTrack()) {
      root_io->CloseTrackFile();
    }
    if (OutputConfig::GetSaveStep()) {
      root_io->CloseStepFile();
    }
  }

  // Print results
  G4cout << "  The run was " << run->GetNumberOfEvent() << " events " << G4endl;

  timer->Stop();
  G4cout << " time:  " << *timer << G4endl;
}
