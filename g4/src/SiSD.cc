#include "SiSD.hh"

#include "DetectorChannel.hh"
#include "SiArrayConfig.hh"
#include "SiDetector.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiSD::SiSD(const G4String& name, const G4String& front_collection_name, const G4String& back_collection_name)
    : G4VSensitiveDetector(name)
{
  collectionName.insert(front_collection_name);
  collectionName.insert(back_collection_name);
  collectionName.insert("SiDepositHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::Initialize(G4HCofThisEvent* hce)
{
  front_hits_collection = new SiHitsCollection(SensitiveDetectorName, collectionName[0]);
  back_hits_collection = new SiHitsCollection(SensitiveDetectorName, collectionName[1]);
  deposit_hits_collection = new SiHitsCollection(SensitiveDetectorName, collectionName[2]);

  auto sd_manager = G4SDManager::GetSDMpointer();
  if (front_hc_id < 0) front_hc_id = sd_manager->GetCollectionID(front_hits_collection);
  if (back_hc_id < 0) back_hc_id = sd_manager->GetCollectionID(back_hits_collection);
  if (deposit_hc_id < 0) deposit_hc_id = sd_manager->GetCollectionID(deposit_hits_collection);

  hce->AddHitsCollection(front_hc_id, front_hits_collection);
  hce->AddHitsCollection(back_hc_id, back_hits_collection);
  hce->AddHitsCollection(deposit_hc_id, deposit_hits_collection);

  front_hit_index.clear();
  back_hit_index.clear();
  deposit_hit_index.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::AccumulateStripHit(SiHitsCollection* collection,
                              std::map<StripKey, G4int>& hit_index,
                              const StripKey& key,
                              G4int detector_id,
                              SiReadoutSide side,
                              G4int strip_id,
                              G4double e_dep,
                              G4Step* step)
{
  SiHit* hit = nullptr;
  const auto found = hit_index.find(key);
  if (found == hit_index.end()) {
    hit = new SiHit();
    collection->insert(hit);
    hit_index[key] = static_cast<G4int>(collection->GetSize()) - 1;

    hit->SetDetectorId(detector_id);
    hit->SetReadoutSide(static_cast<G4int>(side));
    hit->SetStripId(strip_id);
    hit->SetTime(step->GetPreStepPoint()->GetGlobalTime());
  } else {
    hit = (*collection)[found->second];
  }

  hit->AddEdep(e_dep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::AccumulateDepositHit(const DepositKey& key,
                                G4int detector_id,
                                G4double e_dep,
                                const G4ThreeVector& entry_position,
                                const G4ThreeVector& deposit_position,
                                G4Step* step)
{
  SiHit* hit = nullptr;
  const auto found = deposit_hit_index.find(key);

  if (found == deposit_hit_index.end()) {
    hit = new SiHit();
    deposit_hits_collection->insert(hit);
    deposit_hit_index[key] =
        static_cast<G4int>(deposit_hits_collection->GetSize()) - 1;

    const auto track = step->GetTrack();
    hit->SetDetectorId(detector_id);
    hit->SetTrackId(track->GetTrackID());
    hit->SetParentId(track->GetParentID());
    hit->SetPdg(track->GetDefinition()->GetPDGEncoding());
    hit->SetTime(step->GetPreStepPoint()->GetGlobalTime());
    hit->SetEntryPosition(entry_position);
  } else {
    hit = (*deposit_hits_collection)[found->second];
  }

  hit->AddEdepAtPosition(e_dep, deposit_position);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SiSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  const G4double e_dep = step->GetTotalEnergyDeposit();
  if (e_dep <= 0.) return false;

  const auto pre_step_point = step->GetPreStepPoint();
  const auto touchable = pre_step_point->GetTouchable();
  const auto physical = touchable->GetVolume();

  const G4String physical_name = physical->GetName();
  const G4String suffix = "_phy";
  if (physical_name.size() <= suffix.size()) return false;
  const auto suffix_position = physical_name.size() - suffix.size();
  if (physical_name.compare(suffix_position, suffix.size(), suffix) != 0) return false;
  const G4String detector_name = physical_name.substr(0, suffix_position);

  if (SiDetector::map_name_to_ring_id.find(detector_name) == SiDetector::map_name_to_ring_id.end()) return false;

  // The physical-volume copy number is used only as an opaque, globally unique
  // detector id.  Strip ids remain independent electronic-channel indices.
  const G4int detector_id = physical->GetCopyNo();
  const DetectorChannel base = DecodeDetectorCopyNo(detector_id);
  if (base.detector_type != static_cast<G4int>(DetectorType::Si) || base.module_id < 0) return false;

  const G4ThreeVector world_pos = pre_step_point->GetPosition();
  const G4ThreeVector local_pos =
      touchable->GetHistory()->GetTopTransform().TransformPoint(world_pos);

  const G4ThreeVector post_world_pos =
      step->GetPostStepPoint()->GetPosition();
  const G4ThreeVector deposit_position =
      0.5 * (world_pos + post_world_pos);

  const G4int track_id = step->GetTrack()->GetTrackID();
  const DepositKey deposit_key(detector_id, track_id);
  AccumulateDepositHit(
      deposit_key,
      detector_id,
      e_dep,
      world_pos,
      deposit_position,
      step);

  G4int front_strip_id = 0;
  G4int back_strip_id = 0;

  const auto box_it = SiDetector::map_si_box_par.find(detector_name);
  if (box_it != SiDetector::map_si_box_par.end()) {
    // W1: front reads local x; back reads local y.
    const G4double half_x = box_it->second[0] / 2. * mm;
    const G4double half_y = box_it->second[1] / 2. * mm;
    front_strip_id = SiArrayConfig::BarrelPhiStripId(local_pos.x(), half_x);
    back_strip_id = SiArrayConfig::BarrelZStripId(local_pos.y(), half_y);
  } else if (const auto wedge_it = SiDetector::map_si_wedge_par.find(detector_name);
             wedge_it != SiDetector::map_si_wedge_par.end()) {
    // Lampshade wedge: front reads transverse sectors whose physical width
    // grows from the narrow edge to the wide edge; back reads longitudinal
    // rings along local y.
    const G4double inner_half_width = wedge_it->second[0] / 2. * mm;
    const G4double outer_half_width = wedge_it->second[1] / 2. * mm;
    const G4double half_length = wedge_it->second[2] / 2. * mm;
    front_strip_id = SiArrayConfig::LampshadeSectorStripId(
        local_pos.x(),
        local_pos.y(),
        half_length,
        inner_half_width,
        outer_half_width);
    back_strip_id =
        SiArrayConfig::LampshadeRingStripId(local_pos.y(), half_length);
  } else if (const auto tub_it = SiDetector::map_si_par.find(detector_name);
             tub_it != SiDetector::map_si_par.end()) {
    // S3: front reads angular sectors; back reads radial rings.
    const G4double r_outer = tub_it->second[0] / 2. * mm;
    const auto inner_it = SiDetector::map_si_inner_radius.find(detector_name);
    const G4double r_inner = inner_it != SiDetector::map_si_inner_radius.end() ? inner_it->second * mm : 0.;
    front_strip_id = SiArrayConfig::AnnularSectorStripId(local_pos.x(), local_pos.y());
    back_strip_id = SiArrayConfig::AnnularRingStripId(local_pos.x(), local_pos.y(), r_inner, r_outer);
  } else {
    return false;
  }

  const StripKey front_key(detector_id, front_strip_id);
  const StripKey back_key(detector_id, back_strip_id);

  // The same physical deposit induces equal calibrated energy signals on the
  // two DSSD faces.  Each face is accumulated independently.  Contributions
  // from multiple particles sharing one strip merge into a single signal.
  AccumulateStripHit(front_hits_collection, front_hit_index, front_key, detector_id, SiReadoutSide::Front,
                     front_strip_id, e_dep, step);
  AccumulateStripHit(back_hits_collection, back_hit_index, back_key, detector_id, SiReadoutSide::Back,
                     back_strip_id, e_dep, step);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel <= 1) return;

  G4cout << G4endl << "--------> Si front-strip hits: " << front_hits_collection->entries() << G4endl;
  for (G4int i = 0; i < front_hits_collection->entries(); ++i) {
    (*front_hits_collection)[i]->Print();
  }

  G4cout << "--------> Si back-strip hits: " << back_hits_collection->entries() << G4endl;
  for (G4int i = 0; i < back_hits_collection->entries(); ++i) {
    (*back_hits_collection)[i]->Print();
  }

  G4cout << "--------> Si physical hits: "
         << deposit_hits_collection->entries() << G4endl;
  for (G4int i = 0; i < deposit_hits_collection->entries(); ++i) {
    (*deposit_hits_collection)[i]->Print();
  }
}
