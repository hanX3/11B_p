#include "SteppingAction.hh"
#include "RootIO.hh"
#include "OutputConfig.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(RootIO* rio)
    : G4UserSteppingAction(), root_io(rio)
{
  step_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!OutputConfig::GetSaveStep()) return;

  step_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  step_data.track = step->GetTrack()->GetTrackID();
  step_data.de = step->GetTotalEnergyDeposit();
  step_data.pre_x = step->GetPreStepPoint()->GetPosition().x();
  step_data.pre_y = step->GetPreStepPoint()->GetPosition().y();
  step_data.pre_z = step->GetPreStepPoint()->GetPosition().z();
  step_data.pre_total_energy = step->GetPreStepPoint()->GetTotalEnergy();
  step_data.pre_kine_energy = step->GetPreStepPoint()->GetKineticEnergy();
  step_data.post_x = step->GetPostStepPoint()->GetPosition().x();
  step_data.post_y = step->GetPostStepPoint()->GetPosition().y();
  step_data.post_z = step->GetPostStepPoint()->GetPosition().z();
  step_data.post_total_energy = step->GetPostStepPoint()->GetTotalEnergy();
  step_data.post_kine_energy = step->GetPostStepPoint()->GetKineticEnergy();
  step_data.length = step->GetStepLength();
  strcpy(step_data.volume, step->GetTrack()->GetVolume()->GetName());
  strcpy(step_data.particle, step->GetTrack()->GetDefinition()->GetParticleName());
  strcpy(step_data.process, step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());

  /*
  if(strcmp(step_data.particle, "gamma")==0 && step->GetTrack()->GetParentID()==1 && strcmp(step_data.detector, "Target")==0){ root_io->FillStepTree(step_data); }
  */
  root_io->FillStepTree(step_data);
}
