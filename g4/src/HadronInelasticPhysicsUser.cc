#include "HadronInelasticPhysicsUser.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "H11BReaction.hh"
#include "H11BCrossSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HadronInelasticPhysicsUser::HadronInelasticPhysicsUser(const G4String& name)
    : G4VPhysicsConstructor(name) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
HadronInelasticPhysicsUser::~HadronInelasticPhysicsUser() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HadronInelasticPhysicsUser::ConstructParticle() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HadronInelasticPhysicsUser::ConstructProcess()
{
  G4ProcessManager* manager = G4Proton::Proton()->GetProcessManager();

  auto process = new G4HadronInelasticProcess("p_11B", G4Proton::Definition());
  process->RegisterMe(new H11BReaction());
  process->AddDataSet(new H11BCrossSection());
  manager->AddDiscreteProcess(process);
}
