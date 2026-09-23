#include "EventAction.hh"

#include "DetectorChannel.hh"
#include "HPGeSD.hh"
#include "LaBr3SD.hh"
#include "OutputConfig.hh"
#include "RootIO.hh"
#include "SiSD.hh"
#include "VirtualSphereConfig.hh"
#include "VirtualSphereSD.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"

#include <cstddef>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction(RootIO* rio)
    : G4UserEventAction(), hc_id_si_front(-1), hc_id_si_back(-1),
      hc_id_si_deposit(-1), hc_id_hpge(-1), hc_id_labr3(-1),
      hc_id_virtual_sphere(-1), root_io(rio)
{
  event_data.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::~EventAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::BeginOfEventAction(const G4Event*)
{
  const G4bool save_virtual_sphere =
      OutputConfig::GetSaveVirtualSphere() && VirtualSphereConfig::GetEnabled();
  if (!OutputConfig::GetSaveEvent() && !save_virtual_sphere) return;

  auto sd_manager = G4SDManager::GetSDMpointer();
  if (OutputConfig::GetSaveEvent()) {
    if (hc_id_si_front < 0) hc_id_si_front = sd_manager->GetCollectionID("SiSD/SiFrontHitCollection");
    if (hc_id_si_back < 0) hc_id_si_back = sd_manager->GetCollectionID("SiSD/SiBackHitCollection");
    if (hc_id_si_deposit < 0) {
      hc_id_si_deposit =
          sd_manager->GetCollectionID("SiSD/SiDepositHitCollection");
    }
    if (hc_id_hpge < 0) hc_id_hpge = sd_manager->GetCollectionID("HPGeSD/HPGeHitCollection");
    if (hc_id_labr3 < 0) hc_id_labr3 = sd_manager->GetCollectionID("LaBr3SD/LaBr3HitCollection");
  }
  if (save_virtual_sphere && hc_id_virtual_sphere < 0) {
    hc_id_virtual_sphere = sd_manager->GetCollectionID("VirtualSphereSD/VirtualSphereHitCollection");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event* event)
{
  const G4bool save_event = OutputConfig::GetSaveEvent();
  const G4bool save_virtual_sphere = OutputConfig::GetSaveVirtualSphere() && VirtualSphereConfig::GetEnabled();
  if (!save_event && !save_virtual_sphere) return;

  const G4long event_id = event->GetEventID();
  auto hce = event->GetHCofThisEvent();

  if (save_virtual_sphere && hce && hc_id_virtual_sphere >= 0) {
    auto collection = static_cast<VirtualSphereHitsCollection*>(hce->GetHC(hc_id_virtual_sphere));
    if (collection) {
      for (std::size_t i = 0; i < collection->GetSize(); ++i) {
        const auto hit = (*collection)[i];
        if (!hit) continue;

        const auto& momentum = hit->GetMomentum();
        const auto& position = hit->GetPosition();
        const auto direction = momentum.mag() > 0. ? momentum.unit() : G4ThreeVector(0., 0., 1.);
        G4double phi = direction.phi();
        if (phi < 0.) phi += twopi;

        VirtualSphereData sphere_data;
        sphere_data.event_id = event_id;
        sphere_data.track_id = hit->GetTrackId();
        sphere_data.parent_id = hit->GetParentId();
        sphere_data.pdg = hit->GetPdg();
        sphere_data.h11b_reaction_channel = hit->GetH11BReactionChannel();
        sphere_data.h11b_particle_role = hit->GetH11BParticleRole();
        sphere_data.h11b_particle_source = hit->GetH11BParticleSource();
        sphere_data.generator_particle_index = hit->GetGeneratorParticleIndex();
        sphere_data.kinetic_energy_MeV = static_cast<float>(hit->GetKineticEnergy() / MeV);
        sphere_data.px_MeV_c = static_cast<float>(momentum.x() / MeV);
        sphere_data.py_MeV_c = static_cast<float>(momentum.y() / MeV);
        sphere_data.pz_MeV_c = static_cast<float>(momentum.z() / MeV);
        sphere_data.x_mm = static_cast<float>(position.x() / mm);
        sphere_data.y_mm = static_cast<float>(position.y() / mm);
        sphere_data.z_mm = static_cast<float>(position.z() / mm);
        sphere_data.theta_lab_deg = static_cast<float>(direction.theta() / deg);
        sphere_data.phi_lab_deg = static_cast<float>(phi / deg);
        sphere_data.global_time_ns = static_cast<float>(hit->GetGlobalTime() / ns);
        root_io->FillVirtualSphereTree(sphere_data);
      }
    }
  }

  if (save_event) {
    event_data.Clear();
    event_data.event_id = event_id;

    if (hce) {
      const auto append_si_collection = [&](G4int collection_id) {
        if (collection_id < 0) return;
        auto collection = static_cast<SiHitsCollection*>(hce->GetHC(collection_id));
        if (!collection) return;

        for (std::size_t i = 0; i < collection->GetSize(); ++i) {
          const auto hit = (*collection)[i];
          if (!hit || hit->GetEdep() <= 0. || hit->GetDetectorId() < 0 || hit->GetStripId() < 0) continue;

          const DetectorChannel channel = DecodeDetectorCopyNo(hit->GetDetectorId());
          if (channel.detector_type != static_cast<G4int>(DetectorType::Si) || channel.ring_id < 0 ||
              channel.module_id < 0)
            continue;

          event_data.si_subarray_id.push_back(channel.ring_id);
          event_data.si_module_id.push_back(channel.module_id);
          event_data.si_side.push_back(hit->GetReadoutSide());
          event_data.si_strip_id.push_back(hit->GetStripId());
          event_data.si_edep_MeV.push_back(hit->GetEdep() / MeV);
          event_data.si_time_ns.push_back(hit->GetTime() / ns);
        }
      };

      append_si_collection(hc_id_si_front);
      append_si_collection(hc_id_si_back);

      if (hc_id_si_deposit >= 0) {
        auto collection =
            static_cast<SiHitsCollection*>(hce->GetHC(hc_id_si_deposit));

        if (collection) {
          for (std::size_t i = 0; i < collection->GetSize(); ++i) {
            const auto hit = (*collection)[i];
            if (!hit || hit->GetEdep() <= 0. ||
                hit->GetDetectorId() < 0 || hit->GetTrackId() < 0)
              continue;

            const DetectorChannel channel =
                DecodeDetectorCopyNo(hit->GetDetectorId());
            if (channel.detector_type !=
                    static_cast<G4int>(DetectorType::Si) ||
                channel.ring_id < 0 || channel.module_id < 0)
              continue;

            const auto& entry = hit->GetEntryPosition();
            const auto edep_position = hit->GetEdepWeightedPosition();

            event_data.si_hit_detector_id.push_back(hit->GetDetectorId());
            event_data.si_hit_subarray_id.push_back(channel.ring_id);
            event_data.si_hit_module_id.push_back(channel.module_id);
            event_data.si_hit_track_id.push_back(hit->GetTrackId());
            event_data.si_hit_parent_id.push_back(hit->GetParentId());
            event_data.si_hit_pdg.push_back(hit->GetPdg());
            event_data.si_hit_edep_MeV.push_back(hit->GetEdep() / MeV);
            event_data.si_hit_time_ns.push_back(hit->GetTime() / ns);

            event_data.si_hit_x_entry_mm.push_back(entry.x() / mm);
            event_data.si_hit_y_entry_mm.push_back(entry.y() / mm);
            event_data.si_hit_z_entry_mm.push_back(entry.z() / mm);

            event_data.si_hit_x_edep_mm.push_back(
                edep_position.x() / mm);
            event_data.si_hit_y_edep_mm.push_back(
                edep_position.y() / mm);
            event_data.si_hit_z_edep_mm.push_back(
                edep_position.z() / mm);
          }
        }
      }

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

    // Every generated event is written exactly once, including events with no
    // detector signal. Smearing and thresholds remain offline operations.
    root_io->FillEventTree(event_data);
  }

  if (event_id < 10 || event_id % 50000 == 0) {
    G4cout << ">>> Event: " << event_id << G4endl;
  }
}
