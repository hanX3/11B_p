#ifndef SiHit_h
#define SiHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>
#include "tls.hh"

//
class SiHit : public G4VHit
{
public:
  SiHit() = default;
  SiHit(G4int rid, G4int sid);
  SiHit(const SiHit&) = default;
  ~SiHit() override = default;

  // operators
  SiHit& operator=(const SiHit&) = default;
  G4bool operator==(const SiHit&) const;

  inline void *operator new(size_t);
  inline void  operator delete(void*);

  // methods from base class
  void Draw() override;
  void Print() override;

  void SetRingId(G4int rid) { ring_id = rid; }
  void SetSectorId(G4int sid) { sector_id = sid; }
  void SetEdep(G4double de) { e_dep = de; }
  void AddEdep(G4double de) { e_dep += de; }
  void SetPos(G4ThreeVector xyz) { pos = xyz; }

  G4int GetRingId() const { return ring_id; }
  G4int GetSectorId() const { return sector_id; }
  G4double GetEdep() const { return e_dep; }
  G4ThreeVector GetPos() const { return pos; }
  const char* GetDetectorName() const { return "Si"; }

private:
  G4int ring_id = -1; // 1,2
  G4int sector_id = -1; // 
  G4double e_dep;
  G4ThreeVector pos;
};

//
using SiHitsCollection = G4THitsCollection<SiHit>;
extern G4ThreadLocal G4Allocator<SiHit>* SiHitAllocator;

//
inline void* SiHit::operator new(size_t)
{
  if(!SiHitAllocator) SiHitAllocator = new G4Allocator<SiHit>;
  return (void *) SiHitAllocator->MallocSingle();
}

//
inline void SiHit::operator delete(void *hit)
{
  SiHitAllocator->FreeSingle((SiHit*) hit);
}

#endif
