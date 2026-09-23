#ifndef LaBr3Hit_h
#define LaBr3Hit_h 1

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
class LaBr3Hit : public G4VHit
{
public:
  LaBr3Hit() = default;
  LaBr3Hit(G4int rid, G4int sid);
  LaBr3Hit(const LaBr3Hit&) = default;
  ~LaBr3Hit() override = default;

  // operators
  LaBr3Hit& operator=(const LaBr3Hit&) = default;
  G4bool operator==(const LaBr3Hit&) const;

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
  const char* GetDetectorName() const { return "LaBr3"; }

private:
  G4int ring_id = -1; // 1,2
  G4int sector_id = -1; // 
  G4double e_dep;
  G4ThreeVector pos;
};

//
using LaBr3HitsCollection = G4THitsCollection<LaBr3Hit>;
extern G4ThreadLocal G4Allocator<LaBr3Hit>* LaBr3HitAllocator;

//
inline void* LaBr3Hit::operator new(size_t)
{
  if(!LaBr3HitAllocator) LaBr3HitAllocator = new G4Allocator<LaBr3Hit>;
  return (void *) LaBr3HitAllocator->MallocSingle();
}

//
inline void LaBr3Hit::operator delete(void *hit)
{
  LaBr3HitAllocator->FreeSingle((LaBr3Hit*) hit);
}

#endif
