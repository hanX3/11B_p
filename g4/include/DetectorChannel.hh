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

// Decimal copy-number layout (widened to support strip-segmented DSSDs).
//
//   copy_no = type*1e8 + array*1e7 + ring*1e5 + module*1e3 + segment
//
//   field    decimal places   width   max value   used range
//   ------   --------------   -----   ---------   -------------------------
//   type     [1e8]            1       9           1..5 (DetectorType)
//   array    [1e7]            1       9           0
//   ring     [1e5, 1e6]       2       99          1..3
//   module   [1e3, 1e4]       2       99          0..7 (barrel sectors)
//   segment  [1e0..1e2]       3       999         0..255 (16x16 DSSD)
//
// The largest encoded value is type=5 -> 5e8 = 500,000,000, well below the
// 2,147,483,647 G4int (32-bit signed) ceiling.  Because type still occupies
// the most-significant field, `copy_no / kTypeUnit` recovers the detector type
// exactly, and no field can overflow into its neighbour as long as the used
// ranges above are respected.
constexpr G4int kSegmentUnit = 1;         // 3 decimal digits: 0..999
constexpr G4int kModuleUnit = 1000;       // 2 decimal digits: 0..99
constexpr G4int kRingUnit = 100000;       // 2 decimal digits: 0..99
constexpr G4int kArrayUnit = 10000000;    // 1 decimal digit: 0..9
constexpr G4int kTypeUnit = 100000000;    // 1 decimal digit: 1..9
constexpr G4int kMaxSegmentId = 999;      // inclusive upper bound for a segment
constexpr G4int kMaxModuleId = 99;
constexpr G4int kMaxRingId = 99;

inline G4int EncodeDetectorCopyNo(DetectorType detector_type, G4int array_id, G4int ring_id, G4int module_id, G4int segment_id = 0)
{
  return static_cast<G4int>(detector_type) * kTypeUnit + array_id * kArrayUnit + ring_id * kRingUnit +
         module_id * kModuleUnit + segment_id * kSegmentUnit;
}

inline DetectorChannel DecodeDetectorCopyNo(G4int copy_no)
{
  DetectorChannel channel;
  channel.copy_no = copy_no;

  if (copy_no < kTypeUnit) return channel;

  channel.detector_type = copy_no / kTypeUnit;
  G4int remaining = copy_no % kTypeUnit;
  channel.array_id = remaining / kArrayUnit;
  remaining %= kArrayUnit;
  channel.ring_id = remaining / kRingUnit;
  remaining %= kRingUnit;
  channel.module_id = remaining / kModuleUnit;
  channel.segment_id = remaining % kModuleUnit;
  return channel;
}

#endif
