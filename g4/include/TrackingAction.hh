#ifndef TrackingAction_H
#define TrackingAction_H 1

#include "G4UserTrackingAction.hh"

class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Converts secondary-specific H11B creator tags into owned per-track metadata.
// No track ROOT output is produced by this class.
class TrackingAction : public G4UserTrackingAction
{
public:
  TrackingAction() = default;
  ~TrackingAction() override = default;

  void PreUserTrackingAction(const G4Track*) override;
  void PostUserTrackingAction(const G4Track*) override;
};

#endif
