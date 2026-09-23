#include "SiHit.hh"

//
G4ThreadLocal G4Allocator<SiHit>* SiHitAllocator = nullptr;

//
SiHit::SiHit(G4int rid, G4int sid)
: ring_id(rid), sector_id(sid)
{

}

//
G4bool SiHit::operator==(const SiHit& right) const
{
  return (this == &right)?true:false;
}

//
void SiHit::Draw()
{
  G4VVisManager* vis_manager = G4VVisManager::GetConcreteInstance();
  if(vis_manager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    vis_manager->Draw(circle);
  }
}

//
void SiHit::Print()
{
  G4cout
     << " ring id: " << ring_id << " sector_id " << sector_id << " energy dep: " << std::setw(7) << G4BestUnit(e_dep, "Energy")
     << " Position: " << std::setw(7) << G4BestUnit(pos, "Length") << G4endl;
}

