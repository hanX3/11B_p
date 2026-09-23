#ifndef SiHit_H
#define SiHit_H 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4UnitsTable.hh"
#include "G4VHit.hh"

#include <iomanip>
#include "tls.hh"

// Logical DSSD readout sides.  The mapping used by SiSD is:
//   Front: W1 local-x / azimuthal strip; S3 angular-sector strip.
//   Back : W1 local-y / polar strip;     S3 radial-ring strip.
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
  void AddEdep(G4double de)
  {
    e_dep += de;
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

private:
  G4int detector_id = -1;
  G4int readout_side = -1;
  G4int strip_id = -1;
  G4double time_ns = 0.;
  G4double e_dep = 0.;
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
