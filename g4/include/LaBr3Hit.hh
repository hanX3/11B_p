#ifndef LaBr3Hit_H
#define LaBr3Hit_H 1

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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

  inline void* operator new(size_t);
  inline void operator delete(void*);

  // methods from base class
  void Draw() override;
  void Print() override;

  void SetRingId(G4int rid)
  {
    ring_id = rid;
  }
  void SetSectorId(G4int sid)
  {
    sector_id = sid;
  }
  void SetDetectorType(G4int type)
  {
    detector_type = type;
  }
  void SetArrayId(G4int id)
  {
    array_id = id;
  }
  void SetModuleId(G4int id)
  {
    module_id = id;
  }
  void SetSegmentId(G4int id)
  {
    segment_id = id;
  }
  void SetCopyNo(G4int copy)
  {
    copy_no = copy;
  }
  void SetTime(G4double time)
  {
    time_ns = time;
  }
  void SetParticleInfo(G4int pdg_code, G4int track, G4int parent)
  {
    pdg = pdg_code;
    track_id = track;
    parent_id = parent;
  }
  void SetEdep(G4double de)
  {
    e_dep = de;
  }
  void AddEdep(G4double de)
  {
    e_dep += de;
  }
  void SetPos(G4ThreeVector xyz)
  {
    pos = xyz;
  }

  G4int GetRingId() const
  {
    return ring_id;
  }
  G4int GetSectorId() const
  {
    return sector_id;
  }
  G4int GetDetectorType() const
  {
    return detector_type;
  }
  G4int GetArrayId() const
  {
    return array_id;
  }
  G4int GetModuleId() const
  {
    return module_id;
  }
  G4int GetSegmentId() const
  {
    return segment_id;
  }
  G4int GetCopyNo() const
  {
    return copy_no;
  }
  G4double GetTime() const
  {
    return time_ns;
  }
  G4int GetPDG() const
  {
    return pdg;
  }
  G4int GetTrackId() const
  {
    return track_id;
  }
  G4int GetParentId() const
  {
    return parent_id;
  }
  G4double GetEdep() const
  {
    return e_dep;
  }
  G4ThreeVector GetPos() const
  {
    return pos;
  }
  const char* GetDetectorName() const
  {
    return "LaBr3";
  }

private:
  G4int ring_id = -1;   // 1,2
  G4int sector_id = -1; //
  G4int detector_type = 0;
  G4int array_id = 0;
  G4int module_id = -1;
  G4int segment_id = 0;
  G4int copy_no = -1;
  G4int pdg = 0;
  G4int track_id = -1;
  G4int parent_id = -1;
  G4double time_ns = 0.;
  G4double e_dep = 0.;
  G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using LaBr3HitsCollection = G4THitsCollection<LaBr3Hit>;
extern G4ThreadLocal G4Allocator<LaBr3Hit>* LaBr3HitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void* LaBr3Hit::operator new(size_t)
{
  if (!LaBr3HitAllocator)
    LaBr3HitAllocator = new G4Allocator<LaBr3Hit>;
  return (void*)LaBr3HitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void LaBr3Hit::operator delete(void* hit)
{
  LaBr3HitAllocator->FreeSingle((LaBr3Hit*)hit);
}

#endif
