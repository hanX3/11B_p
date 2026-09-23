#ifndef ActionInitialization_H
#define ActionInitialization_H 1

#include "G4VUserActionInitialization.hh"
#include "Rtypes.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ActionInitialization : public G4VUserActionInitialization
{
public:
  explicit ActionInitialization(ULong64_t random_seed);
  ~ActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;

private:
  ULong64_t random_seed = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
