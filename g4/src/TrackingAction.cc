#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "RootIO.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4Positron.hh"


//
TrackingAction::TrackingAction(DetectorConstruction* dc, RootIO* rio)
: G4UserTrackingAction(),
  detector(dc),
  root_io(rio)
{ 

}

//
void TrackingAction::PreUserTrackingAction(const G4Track*)
{ 

}

//
void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  if((MASK&0b010)==0b010){
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
}

