#ifndef ProtonLimiter_h
#define ProtonLimiter_h 1

#include "G4VPhysicsConstructor.hh"

#include "globals.hh"

class ProtonLimiter: public G4VPhysicsConstructor
{
public:
  ProtonLimiter(const G4String& name = "Proton_Limiter");
  virtual ~ProtonLimiter();

protected:
  void ConstructParticle() {};
  void ConstructProcess();
};

#endif
