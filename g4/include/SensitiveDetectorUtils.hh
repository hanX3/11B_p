#ifndef SensitiveDetectorUtils_H
#define SensitiveDetectorUtils_H 1

#include "DetectorChannel.hh"
#include "globals.hh"

#include <map>
#include <optional>

namespace SensitiveDetectorUtils {
struct ChannelAddress
{
  G4int ring_id = -1;
  G4int sector_id = -1;
  G4int hit_index = -1;
  G4int detector_type = 0;
  G4int array_id = 0;
  G4int module_id = -1;
  G4int segment_id = 0;
  G4int copy_no = -1;
};

inline std::optional<ChannelAddress> ResolveChannel(const G4String& physical_name, G4int copy_no, const std::map<G4String, G4int>& ring_ids, const std::map<G4String, G4int>& sector_counts)
{
  const G4String suffix = "_phy";
  if (physical_name.size() <= suffix.size()) return std::nullopt;

  const auto suffix_position = physical_name.size() - suffix.size();
  if (physical_name.compare(suffix_position, suffix.size(), suffix) != 0) return std::nullopt;

  const G4String detector_name = physical_name.substr(0, suffix_position);
  const auto ring_it = ring_ids.find(detector_name);
  const auto sector_it = sector_counts.find(detector_name);
  if (ring_it == ring_ids.end() || sector_it == sector_counts.end()) return std::nullopt;

  const DetectorChannel decoded = DecodeDetectorCopyNo(copy_no);
  G4int sector_id = copy_no;
  G4int detector_type = 0;
  G4int array_id = 0;
  G4int segment_id = 0;
  if (decoded.module_id >= 0) {
    detector_type = decoded.detector_type;
    array_id = decoded.array_id;
    sector_id = decoded.module_id;
    segment_id = decoded.segment_id;
  }

  if (sector_id < 0 || sector_id >= sector_it->second) return std::nullopt;

  G4int hit_index = sector_id;
  for (auto it = sector_counts.begin(); it != sector_it; ++it) {
    hit_index += it->second;
  }

  return ChannelAddress{
      ring_it->second, sector_id, hit_index, detector_type, array_id, sector_id, segment_id, copy_no};
}
} // namespace SensitiveDetectorUtils

#endif
