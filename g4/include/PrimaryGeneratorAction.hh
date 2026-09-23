#ifndef PrimaryGeneratorAction_H
#define PrimaryGeneratorAction_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  G4ParticleGun* GetParticleGun()
  {
    return particle_gun;
  }

  // Set methods
  void SetRandomFlag(G4bool);
  void SetBeamEnergy(G4double b)
  {
    particle_gun->SetParticleEnergy(b);
  }

private:
  G4ParticleGun* particle_gun; // G4 particle gun
};

#endif
