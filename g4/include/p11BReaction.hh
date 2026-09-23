#ifndef p11BReaction_h
#define p11BReaction_h 1

#include "G4HadronicInteraction.hh"
#include "globals.hh"
#include "TH1.h"
class p11BReaction: public G4HadronicInteraction
{
public:
  p11BReaction();
  ~p11BReaction();

public:
  virtual G4HadFinalState *ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &targetNucleus);
  void ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *targetNucleus, G4ParticleDefinition *product0, G4ParticleDefinition *product1, G4ParticleDefinition *product2);
  double_t BreitW(G4double ExAver,G4double T)
  {
      static std::map<std::pair<G4double,G4double>,TH1D*> histCache;
      auto key = std::make_pair(ExAver,T);
      
      if(histCache.find(key) == histCache.end()){
         TH1D* h1 = new TH1D("h1","h1",2*ExAver*100,0,2*ExAver);
         for(int i=0;i<2*ExAver*100;i++)
         {
          G4double x = i/100.;
          G4double F = T*T/4 / ( (x-ExAver)*(x-ExAver) + T*T/4);
          h1->Fill(x,F*1000);
         }
         histCache[key] = h1;
      }
      double_t result = histCache[key]->GetRandom();
      return result;
  }
  void ReactionKinematicGamma(const G4HadProjectile &projectile, G4ParticleDefinition *targetNucleus, G4ParticleDefinition *product0, G4ParticleDefinition *product1);
};


#endif
