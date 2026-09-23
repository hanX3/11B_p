#ifndef H11BReaction_h
#define H11BReaction_h 1

#include "RunAction.hh"
#include "G4String.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"

#include "DataStructure.hh"
#include "RootIO.hh"

//
class H11BReaction: public G4HadronicInteraction
{
public:
  enum ResonanceType{
    Resonance165,
    Resonance675
  };

public:
  H11BReaction(ResonanceType type = Resonance165);
  ~H11BReaction();

public:
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &target);
  void ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *target, G4ParticleDefinition *product1, G4ParticleDefinition *product2, G4ParticleDefinition *product3);

  void ChannelSequential1();
  void ChannelSequential2();
  void ChannelSimultaneous();
  void ChannelGamma();

  void ChooseChannel();

private:
  H11BReactionData reaction_data;
  ResonanceType resonance_type;

private:
  G4double Sample8Be2PlusExcitationEnergy(G4double ex_max);
  G4double Weight8Be2Plus(G4double ex_8Be) const;
  G4bool ChooseAlpha0Channel() const;
  G4double GetAlpha0BranchingRatio() const;
  const char* GetResonanceName() const;
};

#endif
