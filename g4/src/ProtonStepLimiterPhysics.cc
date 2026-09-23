#include "ProtonStepLimiterPhysics.hh"

#include "G4Proton.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

ProtonStepLimiterPhysics::ProtonStepLimiterPhysics(const G4String& name)
    : G4VPhysicsConstructor(name)
{
}

void ProtonStepLimiterPhysics::ConstructParticle()
{
    // Make sure proton is defined.
    G4Proton::ProtonDefinition();
}

void ProtonStepLimiterPhysics::ConstructProcess()
{
    G4ParticleDefinition* proton = G4Proton::ProtonDefinition();

    G4ProcessManager* pmanager = proton->GetProcessManager();
    if (!pmanager) return;

    // Only proton gets the StepLimiter process.
    pmanager->AddDiscreteProcess(new G4StepLimiter());
}
