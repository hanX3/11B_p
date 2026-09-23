#include "SiSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "SiDetector.hh"
#include "SiArrayConfig.hh"
#include "DetectorChannel.hh"
#include <cstring>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiSD::SiSD(const G4String& name, const G4String& hits_collection_name)
    : G4VSensitiveDetector(name)
{
  collectionName.insert(hits_collection_name);

  hits_collection = nullptr;
  hc_id = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  hits_collection = new SiHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  if (hc_id < 0) {
    hc_id = G4SDManager::GetSDMpointer()->GetCollectionID(hits_collection);
  }
  hce->AddHitsCollection(hc_id, hits_collection);

  // Strip-segmented readout: hits are created lazily per fired (module, segment)
  // pair inside ProcessHits, so start each event with an empty lookup table.
  map_copyno_to_hit_index.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool SiSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // energy deposit
  G4double e = step->GetTotalEnergyDeposit();
  if (e == 0.) return false;

  auto pre_step_point = step->GetPreStepPoint();
  auto touchable = pre_step_point->GetTouchable();
  auto physical = touchable->GetVolume();

  // Resolve the detector name from the "<name>_phy" physical volume.
  const G4String physical_name = physical->GetName();
  const G4String suffix = "_phy";
  if (physical_name.size() <= suffix.size()) return false;
  const auto suffix_position = physical_name.size() - suffix.size();
  if (physical_name.compare(suffix_position, suffix.size(), suffix) != 0) return false;
  const G4String detector_name = physical_name.substr(0, suffix_position);

  if (SiDetector::map_name_to_ring_id.find(detector_name) == SiDetector::map_name_to_ring_id.end()) return false;

  // Decode the module-level copy number (segment field is zero on the volume).
  const G4int base_copy_no = physical->GetCopyNo();
  const DetectorChannel base = DecodeDetectorCopyNo(base_copy_no);
  if (base.module_id < 0) return false;

  const G4int detector_type = base.detector_type > 0 ? base.detector_type : static_cast<G4int>(DetectorType::Si);
  const G4int array_id = base.array_id;
  const G4int ring_id = base.ring_id;
  const G4int module_id = base.module_id;

  // Global hit position and its module-local counterpart (method (b): a single
  // sensitive module solid plus local-coordinate strip segmentation).
  const G4ThreeVector world_pos = pre_step_point->GetPosition();
  const G4ThreeVector local_pos = touchable->GetHistory()->GetTopTransform().TransformPoint(world_pos);

  G4int segment_id = 0;
  const auto box_it = SiDetector::map_si_box_par.find(detector_name);
  if (box_it != SiDetector::map_si_box_par.end()) {
    // Barrel module (G4Box): local x across width, local y along the beam axis.
    const G4double half_x = box_it->second[0] / 2. * mm;
    const G4double half_y = box_it->second[1] / 2. * mm;
    segment_id = SiArrayConfig::BarrelSegmentId(local_pos.x(), local_pos.y(), half_x, half_y);
  } else {
    // Annular DSSD (G4Tubs): radial rings + angular sectors in the local x-y plane.
    const auto tub_it = SiDetector::map_si_par.find(detector_name);
    const G4double r_outer = tub_it != SiDetector::map_si_par.end() ? tub_it->second[0] / 2. * mm : 0.;
    const auto inner_it = SiDetector::map_si_inner_radius.find(detector_name);
    const G4double r_inner = inner_it != SiDetector::map_si_inner_radius.end() ? inner_it->second * mm : 0.;
    segment_id = SiArrayConfig::AnnularSegmentId(local_pos.x(), local_pos.y(), r_inner, r_outer);
  }

  // Defensive clamp: the copy-number encoding reserves 3 decimal digits for the
  // segment field, so keep segment ids in range even for extreme strip counts.
  if (segment_id > kMaxSegmentId) segment_id = kMaxSegmentId;
  if (segment_id < 0) segment_id = 0;

  const G4int full_copy_no = EncodeDetectorCopyNo(static_cast<DetectorType>(detector_type), array_id, ring_id, module_id, segment_id);

  // Find-or-create the hit for this (module, segment); each fired strip is an
  // independent hit so that several segments of the same module can coincide.
  SiHit* hit = nullptr;
  auto found = map_copyno_to_hit_index.find(full_copy_no);
  if (found == map_copyno_to_hit_index.end()) {
    hit = new SiHit();
    hits_collection->insert(hit);
    map_copyno_to_hit_index[full_copy_no] = static_cast<G4int>(hits_collection->GetSize()) - 1;

    hit->SetRingId(ring_id);
    hit->SetSectorId(module_id);
    hit->SetDetectorType(detector_type);
    hit->SetArrayId(array_id);
    hit->SetModuleId(module_id);
    hit->SetSegmentId(segment_id);
    hit->SetCopyNo(full_copy_no);

    hit->SetPos(world_pos);
    hit->SetTime(pre_step_point->GetGlobalTime());

    auto track = step->GetTrack();
    hit->SetParticleInfo(track->GetDefinition()->GetPDGEncoding(), track->GetTrackID(), track->GetParentID());
  } else {
    hit = (*hits_collection)[found->second];
  }

  hit->AddEdep(e);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    G4int n_of_hits = hits_collection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << n_of_hits << " hits in the tracker chambers: " << G4endl;
    for (G4int i = 0; i < n_of_hits; i++) {
      (*hits_collection)[i]->Print();
    }
  }
}
