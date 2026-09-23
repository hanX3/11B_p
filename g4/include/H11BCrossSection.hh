#ifndef H11BCrossSection_h
#define H11BCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"

//
class H11BCrossSection: public G4VCrossSectionDataSet
{
public:
  enum ResonanceType{
    Resonance165,
    Resonance675
  };

public:
  H11BCrossSection(ResonanceType type = Resonance165);
  virtual ~H11BCrossSection();
  virtual G4bool IsIsoApplicable(const G4DynamicParticle*, G4int, G4int, const G4Element*, const G4Material*){ return true;}
  virtual G4double GetIsoCrossSection(const G4DynamicParticle *aPar, G4int Z, G4int A, const G4Isotope *iso, const G4Element *elm, const G4Material *mat);

private:
  ResonanceType resonance_type;

private:
  G4double GetEcmValue(G4double projectA, G4double targetA, G4double elab);
  G4double GetSigma165(G4double ecm);
  G4double GetSigma675(G4double ecm);
  G4double GetSigmaBreitWigner(G4double ecm,
                               G4double resonance_energy,
                               G4double total_width,
                               G4double entrance_width_at_resonance,
                               G4double exit_width_at_resonance,
                               G4double spin_stat_factor,
                               G4int entrance_orbital_l);
  G4double GetFullCoulombEntranceWidth(G4double ecm,
                                       G4double resonance_energy,
                                       G4double entrance_width_at_resonance,
                                       G4int entrance_orbital_l);
};

#endif
