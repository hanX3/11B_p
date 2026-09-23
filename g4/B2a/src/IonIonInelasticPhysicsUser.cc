#include "IonIonInelasticPhysicsUser.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "JUNAReaction.hh"
#include "JUNACrossSection.hh"


IonIonInelasticPhysicsUser::IonIonInelasticPhysicsUser(const G4String& name):
  G4VPhysicsConstructor(name)
{
  G4cout << G4endl << "A local inelastic model is activated for all ions" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
IonIonInelasticPhysicsUser::~IonIonInelasticPhysicsUser()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void IonIonInelasticPhysicsUser::ConstructProcess()
{
  JUNAReaction *JUNAmodel = new JUNAReaction();
  JUNACrossSection *JUNAdata = new JUNACrossSection();

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while((*particleIterator)()){
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();

    G4String particleName = particle->GetParticleName();
    if(particleName == "proton"){
      G4HadronInelasticProcess *protonInelasticProcess = new G4HadronInelasticProcess("inelastic",G4Proton::Definition());
      protonInelasticProcess->AddDataSet(JUNAdata);
      protonInelasticProcess->RegisterMe(JUNAmodel);
      pmanager->AddDiscreteProcess(protonInelasticProcess);
    }
  }
}
