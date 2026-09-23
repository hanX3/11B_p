#include "TrackingAction.hh"

#include "RootIO.hh"
#include "OutputConfig.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction(RootIO* rio)
    : G4UserTrackingAction(), root_io(rio) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PreUserTrackingAction(const G4Track*) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  if (!OutputConfig::GetSaveTrack()) return;

  track_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  track_data.track = track->GetTrackID();
  track_data.e = track->GetKineticEnergy();
  track_data.x = track->GetPosition().x();
  track_data.y = track->GetPosition().y();
  track_data.z = track->GetPosition().z();
  track_data.ts = track->GetGlobalTime();
  track_data.length = track->GetTrackLength();
  strcpy(track_data.volume, track->GetVolume()->GetName());
  strcpy(track_data.particle, track->GetDefinition()->GetParticleName());

  root_io->FillTrackTree(track_data);
}
