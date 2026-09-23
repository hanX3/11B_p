#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "PhysListEmStandard.hh"
#include "G4EmStandardPhysics.hh"

#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "IonIonInelasticPhysicsUser.hh"
#include "ProtonLimiter.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysicsList::PhysicsList()
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  RegisterPhysics(new PhysListEmStandard());
  // Default physics
  RegisterPhysics(new G4DecayPhysics(1));
  // Synchroton Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(1));
  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysics(1));
  // Hadron Physics
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT(1));
  RegisterPhysics(new IonIonInelasticPhysicsUser("standard"));
  RegisterPhysics(new ProtonLimiter("standard"));
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(1));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(1));
  // Neutron tracking cut
  RegisterPhysics(new G4NeutronTrackingCut(1));
  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics(1));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
  //SetDefaultCutValue(0.001);
}

