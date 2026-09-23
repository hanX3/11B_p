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

  G4HadronInelasticProcess *process165 = new G4HadronInelasticProcess("p_11B_165", G4Proton::Definition());
  process165->RegisterMe(new H11BReaction(H11BReaction::Resonance165));
  process165->AddDataSet(new H11BCrossSection(H11BCrossSection::Resonance165));
  manager->AddDiscreteProcess(process165);

  G4HadronInelasticProcess *process675 = new G4HadronInelasticProcess("p_11B_675", G4Proton::Definition());
  process675->RegisterMe(new H11BReaction(H11BReaction::Resonance675));
  process675->AddDataSet(new H11BCrossSection(H11BCrossSection::Resonance675));
  manager->AddDiscreteProcess(process675);
}
