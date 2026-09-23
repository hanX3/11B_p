#include "OutputConfig.hh"

#include "G4GenericMessenger.hh"
#include "G4ios.hh"

namespace {
G4bool g_save_reaction = true;
G4bool g_save_event = false;
G4bool g_save_track = false;
G4bool g_save_step = false;
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

G4bool OutputConfig::GetSaveTrack()
{
  return g_save_track;
}

G4bool OutputConfig::GetSaveStep()
{
  return g_save_step;
}

void OutputConfig::SetSaveReaction(G4bool enabled)
{
  g_save_reaction = enabled;
}

void OutputConfig::SetSaveEvent(G4bool enabled)
{
  g_save_event = enabled;
}

void OutputConfig::SetSaveTrack(G4bool enabled)
{
  g_save_track = enabled;
}

void OutputConfig::SetSaveStep(G4bool enabled)
{
  g_save_step = enabled;
}

void OutputConfig::SetSaveReactionCommand(G4bool enabled)
{
  SetSaveReaction(enabled);
}

void OutputConfig::SetSaveEventCommand(G4bool enabled)
{
  SetSaveEvent(enabled);
}

void OutputConfig::SetSaveTrackCommand(G4bool enabled)
{
  SetSaveTrack(enabled);
}

void OutputConfig::SetSaveStepCommand(G4bool enabled)
{
  SetSaveStep(enabled);
}

void OutputConfig::PrintConfigCommand()
{
  G4cout << "\n===== Output runtime configuration =====" << G4endl
         << "  saveReaction = " << (GetSaveReaction() ? "true" : "false") << G4endl
         << "  saveEvent    = " << (GetSaveEvent() ? "true" : "false") << G4endl
         << "  saveTrack    = " << (GetSaveTrack() ? "true" : "false") << G4endl
         << "  saveStep     = " << (GetSaveStep() ? "true" : "false") << G4endl;
}

void OutputConfig::DefineCommands()
{
  messenger = std::make_unique<G4GenericMessenger>(this, "/output/", "ROOT output controls");

  messenger->DeclareMethod("saveReaction", &OutputConfig::SetSaveReactionCommand,
                           "Enable or disable reaction ROOT output");
  messenger->DeclareMethod("saveEvent", &OutputConfig::SetSaveEventCommand,
                           "Enable or disable detector-event ROOT output");
  messenger->DeclareMethod("saveTrack", &OutputConfig::SetSaveTrackCommand,
                           "Enable or disable track ROOT output");
  messenger->DeclareMethod("saveStep", &OutputConfig::SetSaveStepCommand,
                           "Enable or disable step ROOT output");
  messenger->DeclareMethod("printConfig", &OutputConfig::PrintConfigCommand,
                           "Print current ROOT output configuration");
}
