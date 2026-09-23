#ifndef HadronInelasticPhysicsUser_h
#define HadronInelasticPhysicsUser_h 1

#include "G4VPhysicsConstructor.hh"

#include "globals.hh"

//
class HadronInelasticPhysicsUser: public G4VPhysicsConstructor
{
public:
  HadronInelasticPhysicsUser(const G4String &name = "hadron_inelastic_user");
  virtual ~HadronInelasticPhysicsUser();

protected:
  void ConstructParticle() override;
  void ConstructProcess() override;
};

#endif
