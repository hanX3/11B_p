#include "Constants.hh"
#include "DetectorConstruction.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RootIO.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "sys/stat.h" //  mkdir
#include "TROOT.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  gROOT->Reset();

  // mkdir
  mkdir(DATAPATH, 0777);

  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive *ui = 0;
  if(argc==1){
    ui = new G4UIExecutive(argc, argv);
  }

  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);

  //
  G4RunManager *run_manager = new G4RunManager;
  run_manager->SetUserInitialization(new DetectorConstruction());
  run_manager->SetUserInitialization(new PhysicsList());
  run_manager->SetUserAction(new PrimaryGeneratorAction());

  RootIO *root_io = new RootIO();
  run_manager->SetUserAction(new RunAction(root_io));
  run_manager->SetUserAction(new EventAction(new PrimaryGeneratorAction(), root_io));
  if((MASK&0b010)==0b010) run_manager->SetUserAction(new TrackingAction(new DetectorConstruction(), root_io));
  if((MASK&0b100)==0b100) run_manager->SetUserAction(new SteppingAction(new DetectorConstruction(), root_io));
  run_manager->Initialize();

  //
  G4VisManager *vis_manager = new G4VisExecutive;
  vis_manager->Initialize();

  G4UImanager *ui_manager = G4UImanager::GetUIpointer();
  if(!ui){
    G4String command = "/control/execute ";
    G4String file_name = argv[1];
    ui_manager->ApplyCommand(command+file_name);
  }
  else{
    ui_manager->ApplyCommand("/control/execute macros/init_vis.mac");
    if(ui->IsGUI()) {
      ui_manager->ApplyCommand("/control/execute macros/gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  //
  delete vis_manager;
  delete run_manager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
