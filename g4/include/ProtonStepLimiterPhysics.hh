#ifndef ProtonStepLimiterPhysics_h
#define ProtonStepLimiterPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class ProtonStepLimiterPhysics : public G4VPhysicsConstructor
{
public:
    ProtonStepLimiterPhysics(const G4String& name = "ProtonStepLimiterPhysics");
    virtual ~ProtonStepLimiterPhysics() = default;

    virtual void ConstructParticle() override;
    virtual void ConstructProcess() override;
};

#endif
