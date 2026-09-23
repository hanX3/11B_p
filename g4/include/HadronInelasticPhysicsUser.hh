#ifndef HadronInelasticPhysicsUser_H
#define HadronInelasticPhysicsUser_H 1

#include "G4VPhysicsConstructor.hh"

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class HadronInelasticPhysicsUser : public G4VPhysicsConstructor
{
public:
  HadronInelasticPhysicsUser(const G4String& name = "hadron_inelastic_user");
  virtual ~HadronInelasticPhysicsUser();

protected:
  void ConstructParticle() override;
  void ConstructProcess() override;
};

#endif
