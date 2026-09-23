// copy froam extended/radioactivedecay/rdecay02

#ifndef PhysicsList_H
#define PhysicsList_H 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList();
  ~PhysicsList() override = default;

public:
  virtual void SetCuts();
};

#endif
