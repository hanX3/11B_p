#ifndef IonIonInelasticPhysicsUser_h
#define IonIonInelasticPhysicsUser_h 1

#include "G4VPhysicsConstructor.hh"

#include "globals.hh"

class IonIonInelasticPhysicsUser: public G4VPhysicsConstructor
{
public:
  IonIonInelasticPhysicsUser(const G4String& name = "ion_ion_inelastic_user");
  virtual ~IonIonInelasticPhysicsUser();

protected:
  void ConstructParticle() {};
  void ConstructProcess();
};

#endif
