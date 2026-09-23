#include "OutputConfig.hh"

#include "G4GenericMessenger.hh"
#include "G4ios.hh"

namespace {
G4bool g_save_reaction = true;
G4bool g_save_event = false;
G4bool g_save_virtual_sphere = false;
} // namespace

OutputConfig::OutputConfig()
{
  DefineCommands();
}

OutputConfig::~OutputConfig() = default;

G4bool OutputConfig::GetSaveReaction()
{
  return g_save_reaction;
}

G4bool OutputConfig::GetSaveEvent()
{
  return g_save_event;
}

G4bool OutputConfig::GetSaveVirtualSphere()
{
  return g_save_virtual_sphere;
}

void OutputConfig::SetSaveReaction(G4bool enabled)
{
  g_save_reaction = enabled;
}

void OutputConfig::SetSaveEvent(G4bool enabled)
{
  g_save_event = enabled;
}

void OutputConfig::SetSaveVirtualSphere(G4bool enabled)
{
  g_save_virtual_sphere = enabled;
}

void OutputConfig::SetSaveReactionCommand(G4bool enabled)
{
  SetSaveReaction(enabled);
}

void OutputConfig::SetSaveEventCommand(G4bool enabled)
{
  SetSaveEvent(enabled);
}

void OutputConfig::SetSaveVirtualSphereCommand(G4bool enabled)
{
  SetSaveVirtualSphere(enabled);
}

void OutputConfig::PrintConfigCommand()
{
  G4cout << "\n===== Output runtime configuration =====" << G4endl
         << "  saveReaction      = " << (GetSaveReaction() ? "true" : "false") << G4endl
         << "  saveEvent         = " << (GetSaveEvent() ? "true" : "false") << G4endl
         << "  saveVirtualSphere = " << (GetSaveVirtualSphere() ? "true" : "false") << G4endl;
}

void OutputConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/output/", "Unified ROOT output controls");

  messenger->DeclareMethod("saveReaction", &OutputConfig::SetSaveReactionCommand,
                           "Enable or disable the reaction tree");
  messenger->DeclareMethod("saveEvent", &OutputConfig::SetSaveEventCommand,
                           "Enable or disable the detector event tree");
  messenger->DeclareMethod("saveVirtualSphere", &OutputConfig::SetSaveVirtualSphereCommand,
                           "Enable or disable the virtual_sphere tree");
  messenger->DeclareMethod("printConfig", &OutputConfig::PrintConfigCommand,
                           "Print current unified ROOT output configuration");
}
