#ifndef H11BReaction_h
#define H11BReaction_h 1

#include "RunAction.hh"
#include "G4String.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"

#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TSpline.h"

#include "DataStructure.hh"
#include "RootIO.hh"

//
class H11BReaction: public G4HadronicInteraction
{
public:
  H11BReaction();
  ~H11BReaction();

public:
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &target);
  void ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *target, G4ParticleDefinition *product1, G4ParticleDefinition *product2, G4ParticleDefinition *product3);

  void ChannelSequential1();
  void ChannelSequential2();
  void ChannelSimultaneous();
  void ChannelGamma();

  void ChooseChannel();

public:
  static TSpline3 *sp3_seq2over1;

private:
  H11BReactionData reaction_data;

private:
  G4double GetExBreitW(G4double ex_aver, G4double gamma);
  void SetSeq2Over1Spline3();

  G4bool ChooseSeq(G4double e_cm); // 0->seq1 1->seq2
};

#endif
