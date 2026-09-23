#ifndef OutputConfig_H
#define OutputConfig_H 1

#include "globals.hh"

#include <memory>

class G4GenericMessenger;

// Runtime output controls exposed under /output/.
// Defaults are suitable for generator-level reaction validation:
// save reaction ROOT output, skip event/track/step ROOT output.
class OutputConfig
{
public:
  OutputConfig();
  ~OutputConfig();

  static G4bool GetSaveReaction();
  static G4bool GetSaveEvent();
  static G4bool GetSaveTrack();
  static G4bool GetSaveStep();

  static void SetSaveReaction(G4bool enabled);
  static void SetSaveEvent(G4bool enabled);
  static void SetSaveTrack(G4bool enabled);
  static void SetSaveStep(G4bool enabled);

  void SetSaveReactionCommand(G4bool enabled);
  void SetSaveEventCommand(G4bool enabled);
  void SetSaveTrackCommand(G4bool enabled);
  void SetSaveStepCommand(G4bool enabled);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
