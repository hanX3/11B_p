#ifndef DetectorChannel_H
#define DetectorChannel_H 1

#include "globals.hh"

enum class DetectorType
{
  Unknown = 0,
  Si = 1,
  LaBr3 = 2,
  HPGe = 3,
  BGO = 4,
  NaI = 5
};

struct DetectorChannel
{
  G4int detector_type = 0;
  G4int array_id = 0;
  G4int ring_id = -1;
  G4int module_id = -1;
  G4int segment_id = 0;
  G4int copy_no = -1;
};

inline G4int EncodeDetectorCopyNo(DetectorType detector_type, G4int array_id, G4int ring_id, G4int module_id, G4int segment_id = 0)
{
  return static_cast<G4int>(detector_type) * 100000 + array_id * 10000 + ring_id * 1000 + module_id * 10 +
         segment_id;
}

inline DetectorChannel DecodeDetectorCopyNo(G4int copy_no)
{
  DetectorChannel channel;
  channel.copy_no = copy_no;

  if (copy_no < 100000) return channel;

  channel.detector_type = copy_no / 100000;
  G4int remaining = copy_no % 100000;
  channel.array_id = remaining / 10000;
  remaining %= 10000;
  channel.ring_id = remaining / 1000;
  remaining %= 1000;
  channel.module_id = remaining / 10;
  channel.segment_id = remaining % 10;
  return channel;
}

#endif
