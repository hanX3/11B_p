#ifndef JUNACrossSection_h
#define JUNACrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"

class JUNACrossSection: public G4VCrossSectionDataSet
{
public:
  JUNACrossSection();
  virtual ~JUNACrossSection();
  virtual G4bool IsIsoApplicable(const G4DynamicParticle*, G4int, G4int, const G4Element*, const G4Material*){ return true;}
  virtual G4double GetIsoCrossSection(const G4DynamicParticle *aPar, G4int Z, G4int A, const G4Isotope *iso, const G4Element *elm, const G4Material *mat);

private:
  G4double GetEcmValue(G4double projectA, G4double targetA, G4double elab);
  G4double GetSigmaBySfactor(G4double projectZ, G4double projectA, G4double targetZ, G4double targetA, G4double Ecm, G4double Sfactor);
};

#endif
