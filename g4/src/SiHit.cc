#include "SiHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal G4Allocator<SiHit>* SiHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SiHit::operator==(const SiHit& right) const
{
  return (this == &right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiHit::Draw()
{
  // Strip signals have no unique two-dimensional position until a front/back
  // pair has been reconstructed, so they are intentionally not drawn here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiHit::Print()
{
  G4cout << " detector: " << detector_id
         << " side: " << readout_side
         << " strip: " << strip_id
         << " energy dep: " << std::setw(7) << G4BestUnit(e_dep, "Energy")
         << " first time: " << std::setw(7) << G4BestUnit(time_ns, "Time") << G4endl;
}
