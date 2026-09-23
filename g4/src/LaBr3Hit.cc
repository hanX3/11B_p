#include "LaBr3Hit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ThreadLocal G4Allocator<LaBr3Hit>* LaBr3HitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LaBr3Hit::LaBr3Hit(G4int rid, G4int sid)
    : ring_id(rid), sector_id(sid) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool LaBr3Hit::operator==(const LaBr3Hit& right) const
{
  return (this == &right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Hit::Draw()
{
  G4VVisManager* vis_manager = G4VVisManager::GetConcreteInstance();
  if (vis_manager) {
    G4Circle circle(pos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1., 0., 0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    vis_manager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LaBr3Hit::Print()
{
  G4cout << " ring id: " << ring_id << " sector_id " << sector_id << " energy dep: " << std::setw(7) << G4BestUnit(e_dep, "Energy") << " Position: " << std::setw(7) << G4BestUnit(pos, "Length") << G4endl;
}
