#ifndef VirtualSphereHit_H
#define VirtualSphereHit_H 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "tls.hh"

class VirtualSphereHit : public G4VHit
{
public:
  VirtualSphereHit() = default;
  VirtualSphereHit(const VirtualSphereHit&) = default;
  ~VirtualSphereHit() override = default;

  VirtualSphereHit& operator=(const VirtualSphereHit&) = default;
  G4bool operator==(const VirtualSphereHit&) const;

  inline void* operator new(size_t);
  inline void operator delete(void*);

  void Draw() override;
  void Print() override;

  void SetTrackId(G4int value) { track_id = value; }
  void SetParentId(G4int value) { parent_id = value; }
  void SetPdg(G4int value) { pdg = value; }
  void SetH11BReactionChannel(G4int value) { h11b_reaction_channel = value; }
  void SetH11BParticleRole(G4int value) { h11b_particle_role = value; }
  void SetH11BParticleSource(G4int value) { h11b_particle_source = value; }
  void SetGeneratorParticleIndex(G4int value) { generator_particle_index = value; }
  void SetKineticEnergy(G4double value) { kinetic_energy = value; }
  void SetMomentum(const G4ThreeVector& value) { momentum = value; }
  void SetPosition(const G4ThreeVector& value) { position = value; }
  void SetGlobalTime(G4double value) { global_time = value; }

  G4int GetTrackId() const { return track_id; }
  G4int GetParentId() const { return parent_id; }
  G4int GetPdg() const { return pdg; }
  G4int GetH11BReactionChannel() const { return h11b_reaction_channel; }
  G4int GetH11BParticleRole() const { return h11b_particle_role; }
  G4int GetH11BParticleSource() const { return h11b_particle_source; }
  G4int GetGeneratorParticleIndex() const { return generator_particle_index; }
  G4double GetKineticEnergy() const { return kinetic_energy; }
  const G4ThreeVector& GetMomentum() const { return momentum; }
  const G4ThreeVector& GetPosition() const { return position; }
  G4double GetGlobalTime() const { return global_time; }

private:
  G4int track_id = -1;
  G4int parent_id = -1;
  G4int pdg = 0;
  G4int h11b_reaction_channel = 0;
  G4int h11b_particle_role = 0;
  G4int h11b_particle_source = 0;
  G4int generator_particle_index = -1;
  G4double kinetic_energy = 0.;
  G4ThreeVector momentum;
  G4ThreeVector position;
  G4double global_time = 0.;
};

using VirtualSphereHitsCollection = G4THitsCollection<VirtualSphereHit>;
extern G4ThreadLocal G4Allocator<VirtualSphereHit>* VirtualSphereHitAllocator;

inline void* VirtualSphereHit::operator new(size_t)
{
  if (!VirtualSphereHitAllocator) VirtualSphereHitAllocator = new G4Allocator<VirtualSphereHit>;
  return static_cast<void*>(VirtualSphereHitAllocator->MallocSingle());
}

inline void VirtualSphereHit::operator delete(void* hit)
{
  VirtualSphereHitAllocator->FreeSingle(static_cast<VirtualSphereHit*>(hit));
}

#endif
