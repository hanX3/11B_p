#include "EventAction.hh"

#include "HPGeSD.hh"
#include "LaBr3SD.hh"
#include "OutputConfig.hh"
#include "RootIO.hh"
#include "SiSD.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <cstddef>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction(RootIO* rio)
    : G4UserEventAction(), hc_id_si_front(-1), hc_id_si_back(-1), hc_id_hpge(-1), hc_id_labr3(-1), root_io(rio)
{
  event_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{
  if (!OutputConfig::GetSaveEvent()) return;

  auto sd_manager = G4SDManager::GetSDMpointer();
  if (hc_id_si_front < 0) hc_id_si_front = sd_manager->GetCollectionID("SiSD/SiFrontHitCollection");
  if (hc_id_si_back < 0) hc_id_si_back = sd_manager->GetCollectionID("SiSD/SiBackHitCollection");
  if (hc_id_hpge < 0) hc_id_hpge = sd_manager->GetCollectionID("HPGeSD/HPGeHitCollection");
  if (hc_id_labr3 < 0) hc_id_labr3 = sd_manager->GetCollectionID("LaBr3SD/LaBr3HitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* event)
{
  if (!OutputConfig::GetSaveEvent()) return;

  event_data.Clear();
  event_data.event_id = event->GetEventID();

  auto hce = event->GetHCofThisEvent();
  if (hce) {
    const auto append_si_collection = [&](G4int collection_id) {
      if (collection_id < 0) return;
      auto collection = static_cast<SiHitsCollection*>(hce->GetHC(collection_id));
      if (!collection) return;

      for (std::size_t i = 0; i < collection->GetSize(); ++i) {
        const auto hit = (*collection)[i];
        if (!hit || hit->GetEdep() <= 0. || hit->GetDetectorId() < 0 || hit->GetStripId() < 0) continue;

        event_data.si_detector_id.push_back(hit->GetDetectorId());
        event_data.si_side.push_back(hit->GetReadoutSide());
        event_data.si_strip_id.push_back(hit->GetStripId());
        event_data.si_edep_MeV.push_back(hit->GetEdep() / MeV);
        event_data.si_time_ns.push_back(hit->GetTime() / ns);
      }
    };

    append_si_collection(hc_id_si_front);
    append_si_collection(hc_id_si_back);

    if (hc_id_labr3 >= 0) {
      auto collection = static_cast<LaBr3HitsCollection*>(hce->GetHC(hc_id_labr3));
      if (collection) {
        for (std::size_t i = 0; i < collection->GetSize(); ++i) {
          const auto hit = (*collection)[i];
          if (!hit || hit->GetEdep() <= 0. || hit->GetCopyNo() < 0) continue;

          event_data.labr3_detector_id.push_back(hit->GetCopyNo());
          event_data.labr3_edep_MeV.push_back(hit->GetEdep() / MeV);
          event_data.labr3_time_ns.push_back(hit->GetTime() / ns);
        }
      }
    }

    if (hc_id_hpge >= 0) {
      auto collection = static_cast<HPGeHitsCollection*>(hce->GetHC(hc_id_hpge));
      if (collection) {
        for (std::size_t i = 0; i < collection->GetSize(); ++i) {
          const auto hit = (*collection)[i];
          if (!hit || hit->GetEdep() <= 0. || hit->GetCopyNo() < 0) continue;

          event_data.hpge_detector_id.push_back(hit->GetCopyNo());
          event_data.hpge_edep_MeV.push_back(hit->GetEdep() / MeV);
          event_data.hpge_time_ns.push_back(hit->GetTime() / ns);
        }
      }
    }
  }

  // Every generated Geant4 event is written exactly once, including events
  // with no detector signal.  Resolution, thresholds and bad-channel masks are
  // intentionally deferred to offline Python analysis.
  root_io->FillEventTree(event_data);

  const G4int event_id = event->GetEventID();
  if (event_id < 10 || event_id % 50000 == 0) {
    G4cout << ">>> Event: " << event_id << G4endl;
  }
}
