#include "H11BCrossSection.hh"

#include "Constants.hh"
#include "CoulombPenetrability.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <cmath>

//
H11BCrossSection::H11BCrossSection(ResonanceType type)
: resonance_type(type)
{

}

//
H11BCrossSection::~H11BCrossSection()
{

}

//
G4double H11BCrossSection::GetIsoCrossSection(const G4DynamicParticle *aPar, G4int Z, G4int A, const G4Isotope *, const G4Element *, const G4Material *)
{
  G4double crossSection = 0.;
  const G4int aParA = aPar->GetDefinition()->GetBaryonNumber();
  const G4String aParName = aPar->GetDefinition()->GetParticleName();

  if(aParA==1 && aParName=="proton" && Z==5 && A==11){
    G4double ecm = GetEcmValue(aParA, A, aPar->GetKineticEnergy() / keV); //unit:keV

    if(resonance_type == Resonance165){
      crossSection = GetSigma165(ecm);
    }
    else if(resonance_type == Resonance675){
      crossSection = GetSigma675(ecm);
    }

    crossSection *= 1e5;
  }
  else{
    crossSection = 0.;
  }

  crossSection = crossSection*cm2;
  // G4cout<<crossSection<<G4endl;
  return crossSection;
}

//
G4double H11BCrossSection::GetEcmValue(G4double projectA, G4double targetA, G4double elab)
{
  return targetA/(projectA+targetA)*elab;
}

//
G4double H11BCrossSection::GetSigma165(G4double ecm)
{
  if(ecm<1. || ecm>400.) return 0.;

  G4double exit_width = H11B165Alpha0Width/keV + H11B165Alpha1Width/keV;
  if(H11BIncludeGammaChannelInCrossSection){
    exit_width += H11B165GammaWidth/keV;
  }

  return GetSigmaBreitWigner(ecm,
                             H11B165ResonanceEnergy/keV,
                             H11B165TotalWidth/keV,
                             H11B165ProtonWidth/keV,
                             exit_width,
                             H11B165SpinStatFactor,
                             H11B165EntranceOrbitalL);
}

//
G4double H11BCrossSection::GetSigma675(G4double ecm)
{
  if(ecm<1. || ecm>3500.) return 0.;

  G4double exit_width = H11B675Alpha0Width/keV + H11B675Alpha1Width/keV;
  if(H11BIncludeGammaChannelInCrossSection){
    exit_width += H11B675GammaWidth/keV;
  }

  return GetSigmaBreitWigner(ecm,
                             H11B675ResonanceEnergy/keV,
                             H11B675TotalWidth/keV,
                             H11B675ProtonWidth/keV,
                             exit_width,
                             H11B675SpinStatFactor,
                             H11B675EntranceOrbitalL);
}

//
G4double H11BCrossSection::GetSigmaBreitWigner(G4double ecm,
                                               G4double resonance_energy,
                                               G4double total_width,
                                               G4double entrance_width_at_resonance,
                                               G4double exit_width_at_resonance,
                                               G4double spin_stat_factor,
                                               G4int entrance_orbital_l)
{
  if(ecm<=0. || resonance_energy<=0. || total_width<=0. ||
     entrance_width_at_resonance<=0. || exit_width_at_resonance<=0.) return 0.;

  // Single-level Breit-Wigner formula with an energy-dependent p+11B
  // entrance width:
  //
  //   Gamma_p(E) = Gamma_p(Er) * P_L(E,a) / P_L(Er,a)
  //
  // where P_L is the full Coulomb penetrability
  //
  //   P_L = rho / [F_L(eta,rho)^2 + G_L(eta,rho)^2].
  //
  // The ratio is evaluated from precomputed Coulomb-function tables.  The
  // remaining non-entrance width is kept constant, so the total width at Er
  // is unchanged.
  const G4double projectA = 1.;
  const G4double targetA = 11.;
  const G4double amu_c2_keV = 931494.10242;
  const G4double hbarc_keV_fm = 197326.9804;
  const G4double reduced_mass_keV = projectA*targetA/(projectA+targetA)*amu_c2_keV;

  const G4double pi_over_k2_fm2 = pi*hbarc_keV_fm*hbarc_keV_fm/(2.*reduced_mass_keV*ecm);

  G4double entrance_width_e = entrance_width_at_resonance;
  if(H11BUseFullCoulombPenetrabilityInCrossSection){
    entrance_width_e = GetFullCoulombEntranceWidth(ecm,
                                                   resonance_energy,
                                                   entrance_width_at_resonance,
                                                   entrance_orbital_l);
  }

  const G4double missing_width =
      std::max(0., total_width - entrance_width_at_resonance - exit_width_at_resonance);
  const G4double total_width_e = entrance_width_e + exit_width_at_resonance + missing_width;

  const G4double denominator = (ecm-resonance_energy)*(ecm-resonance_energy)
                             + total_width_e*total_width_e/4.;
  const G4double bw_factor = spin_stat_factor*entrance_width_e*exit_width_at_resonance/denominator;

  return pi_over_k2_fm2*bw_factor*1e-26; // fm2 -> cm2
}

//
G4double H11BCrossSection::GetFullCoulombEntranceWidth(G4double ecm,
                                                       G4double resonance_energy,
                                                       G4double entrance_width_at_resonance,
                                                       G4int entrance_orbital_l)
{
  if(ecm<=0. || resonance_energy<=0. || entrance_width_at_resonance<=0.) return 0.;

  const G4double ratio = H11BCoulomb::P11BRatio(ecm,
                                                resonance_energy,
                                                entrance_orbital_l);

  return entrance_width_at_resonance*ratio;
}
