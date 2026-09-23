#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "SiSD.hh"
#include "HPGeSD.hh"
#include "LaBr3SD.hh"
#include "RootIO.hh"
#include "OutputConfig.hh"

#include "TMath.h"

#include <cstddef>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction(RootIO* rio)
    : G4UserEventAction(), root_io(rio)
{
  hc_id_si = -1;
  threshold_si = SiEnergyThreshold;
  energy_resolution_si = SiEnergyResolution;

  hc_id_hpge = -1;
  threshold_hpge = HPGeEnergyThreshold;
  energy_resolution_hpge = HPGeEnergyResolution;

  hc_id_labr3 = -1;
  threshold_labr3 = LaBr3EnergyThreshold;
  energy_resolution_labr3 = LaBr3EnergyResolution;

  event_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{
  if (!OutputConfig::GetSaveEvent()) return;

  auto sd_manager = G4SDManager::GetSDMpointer();
  hc_id_si = sd_manager->GetCollectionID("SiSD/SiHitCollection");
  hc_id_hpge = sd_manager->GetCollectionID("HPGeSD/HPGeHitCollection");
  hc_id_labr3 = sd_manager->GetCollectionID("LaBr3SD/LaBr3HitCollection");

  // G4cout << "Collection IDs cached: " << hc_id_si << " " << hc_id_hpge << " " << hc_id_labr3 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* event)
{
  if (!OutputConfig::GetSaveEvent()) return;

  auto hce = event->GetHCofThisEvent();

    // Si array, ring 1 sector 1
    auto hc_si = static_cast<SiHitsCollection*>(hce->GetHC(hc_id_si));
    for (std::size_t i = 0; i < hc_si->GetSize(); ++i) {
      if ((*hc_si)[i]->GetEdep() <= 0. || (*hc_si)[i]->GetRingId() < 0) continue;

      event_data.event = event->GetEventID();
      event_data.detector_type = (*hc_si)[i]->GetDetectorType();
      event_data.array_id = (*hc_si)[i]->GetArrayId();
      event_data.ring_id = (*hc_si)[i]->GetRingId();
      event_data.module_id = (*hc_si)[i]->GetModuleId();
      event_data.segment_id = (*hc_si)[i]->GetSegmentId();
      event_data.copy_no = (*hc_si)[i]->GetCopyNo();
      event_data.ring = (*hc_si)[i]->GetRingId();
      event_data.sector = (*hc_si)[i]->GetSectorId();
      event_data.e = (*hc_si)[i]->GetEdep();
      event_data.time = (*hc_si)[i]->GetTime();
      event_data.x = (*hc_si)[i]->GetPos().x();
      event_data.y = (*hc_si)[i]->GetPos().y();
      event_data.z = (*hc_si)[i]->GetPos().z();
      event_data.pdg = (*hc_si)[i]->GetPDG();
      event_data.track_id = (*hc_si)[i]->GetTrackId();
      event_data.parent_id = (*hc_si)[i]->GetParentId();
      strcpy(event_data.detector, (*hc_si)[i]->GetDetectorName());

      GausEnergy(energy_resolution_si);
      if (IfThresholdTrigger(threshold_si)) {
        root_io->FillEventTree(event_data);
      }
    }

    // HPGe array, ring 1,2,... sector 1,2,...
    auto hc_hpge = static_cast<HPGeHitsCollection*>(hce->GetHC(hc_id_hpge));
    for (std::size_t i = 0; i < hc_hpge->GetSize(); ++i) {
      if ((*hc_hpge)[i]->GetEdep() <= 0. || (*hc_hpge)[i]->GetRingId() < 0) continue;

      event_data.event = event->GetEventID();
      event_data.detector_type = (*hc_hpge)[i]->GetDetectorType();
      event_data.array_id = (*hc_hpge)[i]->GetArrayId();
      event_data.ring_id = (*hc_hpge)[i]->GetRingId();
      event_data.module_id = (*hc_hpge)[i]->GetModuleId();
      event_data.segment_id = (*hc_hpge)[i]->GetSegmentId();
      event_data.copy_no = (*hc_hpge)[i]->GetCopyNo();
      event_data.ring = (*hc_hpge)[i]->GetRingId();
      event_data.sector = (*hc_hpge)[i]->GetSectorId();
      event_data.e = (*hc_hpge)[i]->GetEdep();
      event_data.time = (*hc_hpge)[i]->GetTime();
      event_data.x = (*hc_hpge)[i]->GetPos().x();
      event_data.y = (*hc_hpge)[i]->GetPos().y();
      event_data.z = (*hc_hpge)[i]->GetPos().z();
      event_data.pdg = (*hc_hpge)[i]->GetPDG();
      event_data.track_id = (*hc_hpge)[i]->GetTrackId();
      event_data.parent_id = (*hc_hpge)[i]->GetParentId();
      strcpy(event_data.detector, (*hc_hpge)[i]->GetDetectorName());

      GausEnergy(energy_resolution_hpge);
      if (IfThresholdTrigger(threshold_hpge)) {
        root_io->FillEventTree(event_data);
      }
    }

    // LaBr3 array, ring 1,2,... sector 1,2,...
    auto hc_labr3 = static_cast<LaBr3HitsCollection*>(hce->GetHC(hc_id_labr3));
    for (std::size_t i = 0; i < hc_labr3->GetSize(); ++i) {
      if ((*hc_labr3)[i]->GetEdep() <= 0. || (*hc_labr3)[i]->GetRingId() < 0) continue;

      event_data.event = event->GetEventID();
      event_data.detector_type = (*hc_labr3)[i]->GetDetectorType();
      event_data.array_id = (*hc_labr3)[i]->GetArrayId();
      event_data.ring_id = (*hc_labr3)[i]->GetRingId();
      event_data.module_id = (*hc_labr3)[i]->GetModuleId();
      event_data.segment_id = (*hc_labr3)[i]->GetSegmentId();
      event_data.copy_no = (*hc_labr3)[i]->GetCopyNo();
      event_data.ring = (*hc_labr3)[i]->GetRingId();
      event_data.sector = (*hc_labr3)[i]->GetSectorId();
      event_data.e = (*hc_labr3)[i]->GetEdep();
      event_data.time = (*hc_labr3)[i]->GetTime();
      event_data.x = (*hc_labr3)[i]->GetPos().x();
      event_data.y = (*hc_labr3)[i]->GetPos().y();
      event_data.z = (*hc_labr3)[i]->GetPos().z();
      event_data.pdg = (*hc_labr3)[i]->GetPDG();
      event_data.track_id = (*hc_labr3)[i]->GetTrackId();
      event_data.parent_id = (*hc_labr3)[i]->GetParentId();
      strcpy(event_data.detector, (*hc_labr3)[i]->GetDetectorName());

      GausEnergy(energy_resolution_labr3);
      if (IfThresholdTrigger(threshold_labr3)) {
        root_io->FillEventTree(event_data);
      }
    }

  // periodic printing
  G4int event_id = event->GetEventID();
  if (event_id < 10 || event_id % 50000 == 0) {
    G4cout << ">>> Event: " << event_id << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::GausEnergy(G4double res)
{
  event_data.e = G4RandGauss::shoot(event_data.e, res * event_data.e / 2.355);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool EventAction::IfThresholdTrigger(G4double threshold)
{
  if (event_data.e >= threshold) return true;
  return false;
}
