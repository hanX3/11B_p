#include "ActionInitialization.hh"
#include "Constants.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "OutputConfig.hh"
#include "OutputPath.hh"
#include "SiArrayConfig.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "Randomize.hh"
#include "Rtypes.h"
#include "TROOT.h"

#include <chrono>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

namespace {
ULong64_t MakeRandomSeed()
{
  const auto now = static_cast<ULong64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

  std::random_device random_device;
  const auto random0 = static_cast<ULong64_t>(random_device());
  const auto random1 = static_cast<ULong64_t>(random_device());

  ULong64_t mixed = now;
  mixed ^= (random0 << 1);
  mixed ^= (random1 << 33);
  mixed ^= (mixed >> 29);
  mixed *= 0x9E3779B97F4A7C15ULL;
  mixed ^= (mixed >> 32);

  return 1ULL + (mixed % 2147483646ULL);
}

G4int ParseThreadCount(const char* value)
{
  const auto threads = std::stoi(value);
  if (threads <= 0) {
    throw std::out_of_range("thread count must be positive");
  }
  return threads;
}

void PrintUsage(const char* program_name)
{
  G4cerr << "Usage:\n"
         << "  " << program_name << "                         # interactive UI mode\n"
         << "  " << program_name << " <macro.mac>              # batch mode, auto seed, 4 threads\n"
         << "  " << program_name << " <macro.mac> <threads>    # batch mode, auto seed, N threads\n";
}

void MakeDataDirectory()
{
  try {
    HBOutputPath::EnsureDataDirectory();
    G4cout << "----> Data directory = " << HBOutputPath::DataDirectory().string() << G4endl;
  } catch (const std::exception& error) {
    G4cerr << "Failed to create data directory: " << error.what() << G4endl;
    throw;
  }
}
} // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  if (argc > 3) {
    PrintUsage(argv[0]);
    return 1;
  }

  MakeDataDirectory();

  // Detect interactive mode if no arguments are given.
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  const ULong64_t random_seed = MakeRandomSeed();
  CLHEP::HepRandom::setTheSeed(static_cast<long>(random_seed));
  G4cout << "\n----> Random seed = " << random_seed << " (auto)" << G4endl;
  ROOT::EnableThreadSafety();

  G4int n_threads = 4;
  if (!ui && argc > 2) {
    try {
      n_threads = ParseThreadCount(argv[2]);
    } catch (const std::exception& error) {
      G4cerr << "Invalid thread count: " << error.what() << G4endl;
      PrintUsage(argv[0]);
      return 2;
    }
  }

  auto run_manager = ui ? G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly)
                        : G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly, n_threads);
  if (!ui) {
    G4cout << "----> Worker threads = " << n_threads << G4endl;
  }

  auto output_config = std::make_unique<OutputConfig>();
  // Si array geometry / strip-segmentation controls under /si/.  Constructed
  // before /run/initialize so geometry-affecting values can be set from a macro
  // ahead of DetectorConstruction, and runs a copy_no encoding self-test.
  auto si_array_config = std::make_unique<SiArrayConfig>();

  auto detector = new DetectorConstruction();
  run_manager->SetUserInitialization(detector);
  run_manager->SetUserInitialization(new PhysicsList());
  run_manager->SetUserInitialization(new ActionInitialization(random_seed));
  run_manager->Initialize();

  auto vis_manager = new G4VisExecutive();
  vis_manager->Initialize();

  auto ui_manager = G4UImanager::GetUIpointer();
  if (!ui) {
    G4String command = "/control/execute ";
    G4String file_name = argv[1];
    ui_manager->ApplyCommand(command + file_name);
  } else {
    ui_manager->ApplyCommand("/control/execute ../macros/init_vis.mac");
    if (ui->IsGUI()) {
      ui_manager->ApplyCommand("/control/execute ../macros/gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  delete vis_manager;
  delete run_manager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
