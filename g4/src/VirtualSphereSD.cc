#include "VirtualSphereSD.hh"

#include "Constants.hh"
#include "H11BTrackInformation.hh"
#include "VirtualSphereConfig.hh"

#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4ios.hh"

#include <cmath>

VirtualSphereSD::VirtualSphereSD(const G4String& name, const G4String& collection_name)
    : G4VSensitiveDetector(name)
{
  collectionName.insert(collection_name);
}

void VirtualSphereSD::Initialize(G4HCofThisEvent* hce)
{
  hits_collection = new VirtualSphereHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (hc_id < 0) hc_id = G4SDManager::GetSDMpointer()->GetCollectionID(hits_collection);
  hce->AddHitsCollection(hc_id, hits_collection);
  recorded_track_ids.clear();
}

G4bool VirtualSphereSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  if (!VirtualSphereConfig::GetEnabled() || !step) return false;

  const auto pre = step->GetPreStepPoint();
  const auto track = step->GetTrack();
  if (!pre || !track || pre->GetStepStatus() != fGeomBoundary) return false;

  const G4int track_id = track->GetTrackID();
  if (recorded_track_ids.find(track_id) != recorded_track_ids.end()) return false;

  const G4ThreeVector relative_position = pre->GetPosition() - G4ThreeVector(0., 0., TargetZPos);
  const G4ThreeVector direction = pre->GetMomentumDirection();
  if (relative_position.dot(direction) <= 0.) return false;

  const G4double kinetic_energy = pre->GetKineticEnergy();
  if (kinetic_energy < VirtualSphereConfig::GetMinKineticEnergy()) return false;

  const auto definition = track->GetDefinition();
  if (!definition) return false;

  if (!VirtualSphereConfig::GetSaveOpticalPhotons() && definition == G4OpticalPhoton::Definition()) return false;

  const G4int pdg = definition->GetPDGEncoding();
  if (!VirtualSphereConfig::GetSaveElectrons() && (pdg == 11 || pdg == -11)) return false;

  auto hit = new VirtualSphereHit();
  hit->SetTrackId(track_id);
  hit->SetParentId(track->GetParentID());
  hit->SetPdg(pdg);
  if (const auto info =
          dynamic_cast<const H11BTrackInformation*>(track->GetUserInformation())) {
    hit->SetH11BReactionChannel(static_cast<G4int>(info->GetReactionChannel()));
    hit->SetH11BParticleRole(static_cast<G4int>(info->GetParticleRole()));
    hit->SetH11BParticleSource(static_cast<G4int>(info->GetParticleSource()));
    hit->SetGeneratorParticleIndex(info->GetGeneratorParticleIndex());
  }
  hit->SetKineticEnergy(kinetic_energy);
  hit->SetMomentum(pre->GetMomentum());
  hit->SetPosition(pre->GetPosition());
  hit->SetGlobalTime(pre->GetGlobalTime());
  hits_collection->insert(hit);

  recorded_track_ids.insert(track_id);
  return true;
}

void VirtualSphereSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel <= 1 || !hits_collection) return;

  G4cout << G4endl << "--------> Virtual-sphere outward crossings: "
         << hits_collection->entries() << G4endl;
  for (G4int i = 0; i < hits_collection->entries(); ++i) {
    (*hits_collection)[i]->Print();
  }
}
