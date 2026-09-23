#include "TrackingAction.hh"

#include "H11BTrackInformation.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  if (!track || track->GetUserInformation()) return;

  H11BReactionChannel channel = H11BReactionChannel::Unknown;
  H11BParticleRole role = H11BParticleRole::Unknown;
  G4int generator_index = -1;

  if (H11BTrackInformation::DecodeCreatorModelTag(
          track->GetCreatorModelID(), channel, role, generator_index)) {
    track->SetUserInformation(new H11BTrackInformation(
        channel,
        role,
        H11BParticleSource::H11BReactionProduct,
        generator_index));
    return;
  }

  const H11BParticleSource source = track->GetParentID() == 0
      ? H11BParticleSource::BeamParticle
      : H11BParticleSource::TransportSecondary;
  track->SetUserInformation(new H11BTrackInformation(
      H11BReactionChannel::Unknown,
      H11BParticleRole::Unknown,
      source,
      -1));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track*)
{
}
