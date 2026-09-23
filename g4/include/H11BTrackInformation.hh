#ifndef H11BTrackInformation_H
#define H11BTrackInformation_H 1

#include "H11BParticleLabel.hh"

#include "G4VUserTrackInformation.hh"

// Per-track generator metadata. Direct reaction products receive this through
// a secondary-specific creator-model tag, which TrackingAction converts into
// user information before transport begins. Transport daughters get a fresh
// TransportSecondary label and never inherit a generator particle role.
class H11BTrackInformation : public G4VUserTrackInformation
{
public:
  H11BTrackInformation() = default;
  H11BTrackInformation(H11BReactionChannel channel,
                       H11BParticleRole role,
                       H11BParticleSource source,
                       G4int generator_index);
  ~H11BTrackInformation() override = default;

  void Print() const override;

  H11BReactionChannel GetReactionChannel() const { return reaction_channel; }
  H11BParticleRole GetParticleRole() const { return particle_role; }
  H11BParticleSource GetParticleSource() const { return particle_source; }
  G4int GetGeneratorParticleIndex() const { return generator_particle_index; }

  static G4int CreatorModelTag(H11BReactionChannel channel,
                               H11BParticleRole role,
                               G4int generator_index);
  static G4bool DecodeCreatorModelTag(G4int creator_model_id,
                                      H11BReactionChannel& channel,
                                      H11BParticleRole& role,
                                      G4int& generator_index);

private:
  H11BReactionChannel reaction_channel = H11BReactionChannel::Unknown;
  H11BParticleRole particle_role = H11BParticleRole::Unknown;
  H11BParticleSource particle_source = H11BParticleSource::Unknown;
  G4int generator_particle_index = -1;
};

#endif
