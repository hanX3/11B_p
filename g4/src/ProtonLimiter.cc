#include "ProtonLimiter.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"



ProtonLimiter::ProtonLimiter(const G4String& name):
  G4VPhysicsConstructor(name)
{
  G4cout << G4endl << "A local inelastic model is activated for all ions" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ProtonLimiter::~ProtonLimiter()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ProtonLimiter::ConstructProcess()
{

auto particleIterator = GetParticleIterator();
particleIterator->reset();
while((*particleIterator)()){
  G4ParticleDefinition* particle = particleIterator->value();
  G4String name = particle->GetParticleName();
  if(name == "proton"){
    G4ProcessManager* pm = particle->GetProcessManager();
    pm->AddDiscreteProcess(new G4StepLimiter());
  }
}
}
