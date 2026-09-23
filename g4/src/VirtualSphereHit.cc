#include "VirtualSphereHit.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<VirtualSphereHit>* VirtualSphereHitAllocator = nullptr;

G4bool VirtualSphereHit::operator==(const VirtualSphereHit& right) const
{
  return this == &right;
}

void VirtualSphereHit::Draw() {}

void VirtualSphereHit::Print()
{
  G4cout << " virtual-sphere track=" << track_id
         << " parent=" << parent_id
         << " pdg=" << pdg
         << " H11B channel=" << h11b_reaction_channel
         << " role=" << h11b_particle_role
         << " source=" << h11b_particle_source
         << " generator index=" << generator_particle_index
         << " kinetic=" << G4BestUnit(kinetic_energy, "Energy")
         << " position=(" << position.x() << ", " << position.y() << ", " << position.z() << ")"
         << G4endl;
}
