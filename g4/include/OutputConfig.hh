#ifndef OutputConfig_H
#define OutputConfig_H 1

#include "globals.hh"

#include <memory>

class G4GenericMessenger;

// Runtime controls for trees in the unified ROOT output file.
class OutputConfig
{
public:
  OutputConfig();
  ~OutputConfig();

  static G4bool GetSaveReaction();
  static G4bool GetSaveEvent();
  static G4bool GetSaveVirtualSphere();

  static void SetSaveReaction(G4bool enabled);
  static void SetSaveEvent(G4bool enabled);
  static void SetSaveVirtualSphere(G4bool enabled);

  void SetSaveReactionCommand(G4bool enabled);
  void SetSaveEventCommand(G4bool enabled);
  void SetSaveVirtualSphereCommand(G4bool enabled);
  void PrintConfigCommand();

private:
  void DefineCommands();

  std::unique_ptr<G4GenericMessenger> messenger;
};

#endif
