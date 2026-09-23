#include "HPGeSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "HPGeDetector.hh"
#include "SensitiveDetectorUtils.hh"
#include <cstring>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HPGeSD::HPGeSD(const G4String& name, const G4String& hits_collection_name)
    : G4VSensitiveDetector(name)
{
  collectionName.insert(hits_collection_name);

  hits_collection = nullptr;
  hc_id = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  hits_collection = new HPGeHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  if (hc_id < 0) {
    hc_id = G4SDManager::GetSDMpointer()->GetCollectionID(hits_collection);
  }
  hce->AddHitsCollection(hc_id, hits_collection);

  for (auto it = HPGeDetector::map_name_to_sectors.begin(); it != HPGeDetector::map_name_to_sectors.end(); it++) {
    for (auto j = 0; j < it->second; j++) {
      hits_collection->insert(new HPGeHit());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool HPGeSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // energy deposit
  G4double e = step->GetTotalEnergyDeposit();
  if (e == 0.) return false;

  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto physical = touchable->GetVolume();
  auto copy_no = physical->GetCopyNo();
  const auto channel = SensitiveDetectorUtils::ResolveChannel(physical->GetName(), copy_no, HPGeDetector::map_name_to_ring_id, HPGeDetector::map_name_to_sectors);
  if (!channel || channel->hit_index >= static_cast<G4int>(hits_collection->GetSize())) return false;

  /*
  G4cout << "-----> physical name " << det_name << G4endl;
  G4cout << "-----> hc_id " << hc_id << G4endl;
  G4cout << "-----> in HPGeSD ProcessHits function copy_no " << copy_no << G4endl;
  G4cout << "-----> in HPGeSD ProcessHits function ring_id " << ring_id << G4endl;
  G4cout << "-----> in HPGeSD ProcessHits function sector_id " << sector_id << G4endl;
  */

  // check if the first touch
  auto hit = (*hits_collection)[channel->hit_index];
  if (hit->GetRingId() < 0 || hit->GetSectorId() < 0) {
    hit->SetRingId(channel->ring_id);
    hit->SetSectorId(channel->sector_id);
    hit->SetDetectorType(channel->detector_type > 0 ? channel->detector_type : static_cast<G4int>(DetectorType::HPGe));
    hit->SetArrayId(channel->array_id);
    hit->SetModuleId(channel->module_id);
    hit->SetSegmentId(channel->segment_id);
    hit->SetCopyNo(channel->copy_no);

    auto pre_step_point = step->GetPreStepPoint();
    hit->SetPos(pre_step_point->GetPosition());
    hit->SetTime(pre_step_point->GetGlobalTime());

    auto track = step->GetTrack();
    hit->SetParticleInfo(track->GetDefinition()->GetPDGEncoding(), track->GetTrackID(), track->GetParentID());
  }
  hit->AddEdep(e);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HPGeSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    G4int n_of_hits = hits_collection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << n_of_hits << " hits in the tracker chambers: " << G4endl;
    for (G4int i = 0; i < n_of_hits; i++) {
      (*hits_collection)[i]->Print();
    }
  }
}
