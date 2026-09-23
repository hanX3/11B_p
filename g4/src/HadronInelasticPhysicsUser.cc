#include "HadronInelasticPhysicsUser.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "H11BReaction.hh"
#include "H11BCrossSection.hh"

//
HadronInelasticPhysicsUser::HadronInelasticPhysicsUser(const G4String &name)
: G4VPhysicsConstructor(name)
{

}

//
HadronInelasticPhysicsUser::~HadronInelasticPhysicsUser()
{

}

//
void HadronInelasticPhysicsUser::ConstructParticle()
{

}

//
void HadronInelasticPhysicsUser::ConstructProcess()
{
  G4ProcessManager *manager = G4Proton::Proton()->GetProcessManager();
  G4HadronInelasticProcess *process = new G4HadronInelasticProcess("p_11B", G4Proton::Definition());
  
  process->RegisterMe(new H11BReaction());
  process->AddDataSet(new H11BCrossSection());
  manager->AddDiscreteProcess(process);
}

