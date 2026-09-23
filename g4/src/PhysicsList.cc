#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"

#include "HadronInelasticPhysicsUser.hh"

#include "ProtonStepLimiterPhysics.hh"

//
PhysicsList::PhysicsList()
{
  SetVerboseLevel(1);

  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new HadronInelasticPhysicsUser());

  RegisterPhysics(new ProtonStepLimiterPhysics());
}

//
void PhysicsList::SetCuts()
{

}

