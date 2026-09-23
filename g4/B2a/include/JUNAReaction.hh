#ifndef JUNAReaction_h
#define JUNAReaction_h 1

#include "G4HadronicInteraction.hh"
#include "globals.hh"
class JUNAReaction: public G4HadronicInteraction
{
public:
  JUNAReaction();
  ~JUNAReaction();

public:
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &targetNucleus);
  void ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *targetNucleus, G4ParticleDefinition *product0, G4ParticleDefinition *product1, G4ParticleDefinition *product2);
  void ReactionKinematicGamma(const G4HadProjectile &projectile, G4ParticleDefinition *targetNucleus, G4ParticleDefinition *product0, G4ParticleDefinition *product1);
};


#endif
