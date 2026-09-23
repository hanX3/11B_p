#include "PrimaryGeneratorAction.hh"
#include "Constants.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4DynamicParticle.hh"

#include "Randomize.hh"
#include "TMatrixD.h"

//
PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
  G4int n_of_particles = 1;
  particle_gun  = new G4ParticleGun(n_of_particles);

  G4ParticleDefinition *particle_definition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  particle_gun->SetParticleDefinition(particle_definition);

  SetBeamEnergy(BeamEnergy);
}

//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particle_gun;
}

//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* an_event)
{
  G4double r0 = BeamR*G4UniformRand();
  G4double theta = (2.*CLHEP::pi)*G4UniformRand();
  G4double x0 = r0*std::sin(theta);
  G4double y0 = r0*std::cos(theta);
  G4double z0 = 0.;

  particle_gun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  particle_gun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  //set energy
  particle_gun->SetParticleEnergy(particle_energy);

  particle_gun->GeneratePrimaryVertex(an_event);
}