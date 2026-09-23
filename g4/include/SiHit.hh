#ifndef SiHit_H
#define SiHit_H 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4VHit.hh"

#include <iomanip>
#include "tls.hh"

// Logical DSSD readout sides.  The mapping used by SiSD is:
//   Front: W1 local-x; Lampshade transverse sector; S3 angular sector.
//   Back : W1 local-y; Lampshade longitudinal ring; S3 radial ring.
enum class SiReadoutSide : G4int
{
  Front = 0,
  Back = 1
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class SiHit : public G4VHit
{
public:
  SiHit() = default;
  SiHit(const SiHit&) = default;
  ~SiHit() override = default;

  SiHit& operator=(const SiHit&) = default;
  G4bool operator==(const SiHit&) const;

  inline void* operator new(size_t);
  inline void operator delete(void*);

  void Draw() override;
  void Print() override;

  void SetDetectorId(G4int id)
  {
    detector_id = id;
  }
  void SetReadoutSide(G4int side)
  {
    readout_side = side;
  }
  void SetStripId(G4int id)
  {
    strip_id = id;
  }
  void SetTime(G4double time)
  {
    time_ns = time;
  }
  void SetTrackId(G4int id)
  {
    track_id = id;
  }
  void SetParentId(G4int id)
  {
    parent_id = id;
  }
  void SetPdg(G4int value)
  {
    pdg = value;
  }
  void SetEntryPosition(const G4ThreeVector& position)
  {
    if (has_entry_position) return;
    entry_position = position;
    has_entry_position = true;
  }
  void AddEdep(G4double de)
  {
    e_dep += de;
  }
  void AddEdepAtPosition(G4double de, const G4ThreeVector& position)
  {
    if (de <= 0.) return;
    e_dep += de;
    edep_weighted_position += position * de;
  }

  G4int GetDetectorId() const
  {
    return detector_id;
  }
  G4int GetReadoutSide() const
  {
    return readout_side;
  }
  G4int GetStripId() const
  {
    return strip_id;
  }
  G4double GetTime() const
  {
    return time_ns;
  }
  G4double GetEdep() const
  {
    return e_dep;
  }
  G4int GetTrackId() const
  {
    return track_id;
  }
  G4int GetParentId() const
  {
    return parent_id;
  }
  G4int GetPdg() const
  {
    return pdg;
  }
  const G4ThreeVector& GetEntryPosition() const
  {
    return entry_position;
  }
  G4ThreeVector GetEdepWeightedPosition() const
  {
    return e_dep > 0. ? edep_weighted_position / e_dep : G4ThreeVector();
  }

private:
  G4int detector_id = -1;
  G4int readout_side = -1;
  G4int strip_id = -1;
  G4int track_id = -1;
  G4int parent_id = -1;
  G4int pdg = 0;
  G4double time_ns = 0.;
  G4double e_dep = 0.;
  G4bool has_entry_position = false;
  G4ThreeVector entry_position;
  G4ThreeVector edep_weighted_position;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using SiHitsCollection = G4THitsCollection<SiHit>;
extern G4ThreadLocal G4Allocator<SiHit>* SiHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void* SiHit::operator new(size_t)
{
  if (!SiHitAllocator)
    SiHitAllocator = new G4Allocator<SiHit>;
  return (void*)SiHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void SiHit::operator delete(void* hit)
{
  SiHitAllocator->FreeSingle((SiHit*)hit);
}

#endif
