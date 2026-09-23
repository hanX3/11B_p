#include "PrimaryGeneratorAction.hh"
#include "BeamConfig.hh"
#include "Constants.hh"
#include "SiArrayConfig.hh"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  G4int n_of_particles = 1;
  particle_gun = new G4ParticleGun(n_of_particles);

  G4ParticleDefinition* particle_definition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  particle_gun->SetParticleDefinition(particle_definition);

  SetBeamEnergy(BeamEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particle_gun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* an_event)
{
  // Deterministic detector-commissioning mode (default off): emit one alpha
  // from the target centre at a macro-controlled lab direction.  This provides
  // an unambiguous input for validating DSSD strip IDs, including backward S3.
  if (SiArrayConfig::GetFixedAlphaTestBeam()) {
    GenerateFixedAlphaTestPrimary(an_event);
    return;
  }

  // Forward-cone commissioning mode retained for broad forward-S3 illumination.
  if (SiArrayConfig::GetForwardAlphaTestBeam()) {
    GenerateForwardAlphaTestPrimary(an_event);
    return;
  }

  G4double x0 = BeamConfig::GetOffsetX();
  G4double y0 = BeamConfig::GetOffsetY();

  switch (BeamConfig::GetProfile()) {
  case BeamProfile::Point:
    break;

  case BeamProfile::UniformDisk: {
    const G4double radius =
        BeamConfig::GetRadius() * std::sqrt(G4UniformRand());
    const G4double phi = CLHEP::twopi * G4UniformRand();
    x0 += radius * std::cos(phi);
    y0 += radius * std::sin(phi);
    break;
  }

  case BeamProfile::Gaussian:
    x0 += G4RandGauss::shoot(0., BeamConfig::GetSigmaX());
    y0 += G4RandGauss::shoot(0., BeamConfig::GetSigmaY());
    break;
  }

  particle_gun->SetParticlePosition(
      G4ThreeVector(x0, y0, BeamConfig::GetZ()));
  particle_gun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  particle_gun->GeneratePrimaryVertex(an_event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GenerateFixedAlphaTestPrimary(G4Event* an_event)
{
  const G4double theta = SiArrayConfig::GetFixedAlphaTestBeamTheta();
  const G4double phi = SiArrayConfig::GetFixedAlphaTestBeamPhi();
  const G4double sin_theta = std::sin(theta);
  const G4ThreeVector direction(sin_theta * std::cos(phi),
                                sin_theta * std::sin(phi),
                                std::cos(theta));

  particle_gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("alpha"));
  particle_gun->SetParticleEnergy(SiArrayConfig::GetFixedAlphaTestBeamEnergy());
  particle_gun->SetParticlePosition(G4ThreeVector(0., 0., TargetZPos));
  particle_gun->SetParticleMomentumDirection(direction);
  particle_gun->GeneratePrimaryVertex(an_event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GenerateForwardAlphaTestPrimary(G4Event* an_event)
{
  // Launch point just downstream of the target backing (target back face at
  // TargetZPos + TargetThickness/2, backing extends TargetBackingThickness
  // further) so the alpha starts in vacuum and flies toward the forward ring.
  const G4double z0 = TargetZPos + TargetThickness / 2. + TargetBackingThickness + 0.5 * mm;
  const G4double z_ring = TargetZPos + SiArrayConfig::GetForwardDistance();
  const G4double dz = std::max(z_ring - z0, 1. * mm);

  // Sample a forward cone whose half-angles map onto the annular active radius
  // (inner ~8 mm, outer ~60 mm), staying safely inside both edges.
  const G4double r_inner_target = 12. * mm;
  const G4double r_outer_target = 55. * mm;
  const G4double cos_min = dz / std::sqrt(dz * dz + r_outer_target * r_outer_target);
  const G4double cos_max = dz / std::sqrt(dz * dz + r_inner_target * r_inner_target);
  const G4double cos_theta = cos_min + (cos_max - cos_min) * G4UniformRand();
  const G4double sin_theta = std::sqrt(std::max(0., 1. - cos_theta * cos_theta));
  const G4double phi = CLHEP::twopi * G4UniformRand();

  const G4ThreeVector direction(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);

  particle_gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("alpha"));
  particle_gun->SetParticleEnergy(SiArrayConfig::GetForwardAlphaTestBeamEnergy());
  particle_gun->SetParticlePosition(G4ThreeVector(0., 0., z0));
  particle_gun->SetParticleMomentumDirection(direction);
  particle_gun->GeneratePrimaryVertex(an_event);
}
