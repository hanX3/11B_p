#include "H11BReaction.hh"
#include "H11BAngularDistribution.hh"

#include "globals.hh"
#include "G4RunManager.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleMomentum.hh"
#include "G4IonTable.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4AnalysisManager.hh"

#include "TH1.h"
#include "TString.h"
#include <TGenPhaseSpace.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

#include "Constants.hh"
#include "CoulombPenetrability.hh"

#include "Randomize.hh"

//
H11BReaction::H11BReaction(ResonanceType type)
: G4HadronicInteraction(),
  resonance_type(type)
{
  SetMinEnergy(0. * CLHEP::keV);
  SetMaxEnergy(100. * CLHEP::MeV);
  isBlocked = false;
}

//
H11BReaction::~H11BReaction()
{

}

//
G4HadFinalState *H11BReaction::ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &target)
{
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(G4HadFinalStateStatus::stopAndKill);

  G4int z_target = target.GetZ_asInt();
  const G4int a_target = target.GetA_asInt();
  const G4int a_project = projectile.GetDefinition()->GetBaryonNumber();

  if(a_project==1 && a_target==11){

    G4ParticleDefinition *target = G4IonTable::GetIonTable()->GetIon(z_target, a_target, 0. * CLHEP::keV);
    G4ParticleDefinition *product1 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition *product2 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition *product3 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);

    ReactionKinematic(projectile, target, product1, product2, product3);
  }

  return &theParticleChange;
}

//
void H11BReaction::ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *target, G4ParticleDefinition *product1, G4ParticleDefinition *product2, G4ParticleDefinition *product3)
{
  G4double m_projectile = projectile.GetDefinition()->GetPDGMass();
  G4double m_target = target->GetPDGMass();
  
  G4ParticleDefinition *particle_4He = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
  G4double m_4He = particle_4He->GetPDGMass();

  G4ParticleDefinition *particle_8Be = G4IonTable::GetIonTable()->GetIon(4, 8, 0. * CLHEP::keV);
  G4double m_8Be = particle_8Be->GetPDGMass();

  G4double kinetic_lab_projectile = projectile.GetKineticEnergy(); // MeV
  G4double s = m_projectile*m_projectile + m_target*m_target + 2*m_target*(m_projectile + kinetic_lab_projectile);

  //
  // p+11B -> alpha+8Be
  G4double ex_8Be = 0.;
  G4bool flag_seq1 = 0, flag_seq2 = 0, flag_simu = 0;
  if(ChooseAlpha0Channel()) flag_seq2 = 1;
  else flag_seq1 = 1;

  if(flag_seq1){
    G4double ex_max = std::max(0., sqrt(s) - m_4He - m_8Be);
    ex_8Be = Sample8Be2PlusExcitationEnergy(ex_max);
  }

  m_8Be += ex_8Be;
  
  
  G4double e_cm_projectile = (s + m_projectile*m_projectile - m_target*m_target) / (2*sqrt(s));
  G4double e_cm_target = (s + m_target*m_target - m_projectile*m_projectile) / (2*sqrt(s));
  G4double e_cm_alpha1 = (s + m_4He*m_4He - m_8Be*m_8Be) / (2*sqrt(s));
  G4double e_cm_8Be = (s + m_8Be*m_8Be - m_4He*m_4He) / (2*sqrt(s));
  G4double p_cm_projectile = sqrt(e_cm_projectile*e_cm_projectile - m_projectile*m_projectile);
  G4double p_cm_target = sqrt(e_cm_target*e_cm_target - m_target*m_target);
  G4double p_cm_alpha1 = sqrt(e_cm_alpha1*e_cm_alpha1 - m_4He*m_4He);
  G4double p_cm_8Be = sqrt(e_cm_8Be*e_cm_8Be - m_8Be*m_8Be);
  
  G4LorentzVector lv_cm_total(0, 0, 0, sqrt(s));
  
  G4ThreeVector momentum_lab_projectile = projectile.GetMomentumDirection();
  momentum_lab_projectile *= sqrt(kinetic_lab_projectile*kinetic_lab_projectile + 2*kinetic_lab_projectile*m_projectile);
  G4ThreeVector momentum_lab_target(0,0,0);
  
  G4LorentzVector lv_lab_projectile(momentum_lab_projectile, m_projectile+kinetic_lab_projectile);
  G4LorentzVector lv_lab_target(momentum_lab_target, m_target);
  G4LorentzVector lv_lab_total = lv_lab_projectile+lv_lab_target;
  
  //
  // First breakup: 12C* -> primary alpha + 8Be.
  //
  // For the 165-keV 2+ resonance, optionally sample the primary alpha
  // direction relative to the beam axis with a branch-dependent Legendre
  // distribution.  For the 675-keV 2- resonance, keep the previous isotropic
  // primary-alpha direction by default; its important implemented structure is
  // the L=1/L=3 alpha1 angular correlation in the following 8Be(2+) decay.
  G4ThreeVector dir_cm_alpha1;
  if(resonance_type == Resonance165 && H11B165UsePrimaryAngularDistribution){
    dir_cm_alpha1 =
      H11BAngularDistribution::Sample165PrimaryAlphaDirection(projectile.GetMomentumDirection(),
                                                              flag_seq1);
  }
  else{
    dir_cm_alpha1 = H11BAngularDistribution::SampleIsotropicDirection();
  }

  G4double theta_cm_alpha1 = dir_cm_alpha1.theta();
  G4double phi_cm_alpha1 = dir_cm_alpha1.phi();

  G4LorentzVector lv_cm_alpha1(p_cm_alpha1*dir_cm_alpha1.x(),
                               p_cm_alpha1*dir_cm_alpha1.y(),
                               p_cm_alpha1*dir_cm_alpha1.z(),
                               e_cm_alpha1);
  G4LorentzVector lv_cm_8Be = lv_cm_total - lv_cm_alpha1;
  
  G4LorentzVector lv_lab_alpha1 = lv_cm_alpha1;
  lv_lab_alpha1.boost(lv_lab_total.boostVector());

  G4LorentzVector lv_lab_8Be = lv_cm_8Be;
  lv_lab_8Be.boost(lv_lab_total.boostVector());

#ifdef ReactionCout 
  G4cout << "mass projectile " << m_projectile << G4endl;
  G4cout << "mass target " << m_target << G4endl;
  G4cout << "mass 4He " << m_4He << G4endl;
  G4cout << "mass 8Be " << m_8Be << G4endl;
  G4cout << "kinetic_lab_projectile : " << kinetic_lab_projectile << G4endl;
  G4cout << "ex_8Be : " << ex_8Be << " MeV" << G4endl;
  G4cout << "kinetic_energy_cm_projectile : " << e_cm_projectile-m_projectile << " MeV" << G4endl;  
  G4cout << "kinetic_energy_cm_target : " << e_cm_target-m_target << " MeV" << G4endl;  
  G4cout << "kinetic_energy_cm_alpha1 : " << e_cm_alpha1-m_4He << " MeV" << G4endl;  
  G4cout << "kinetic_energy_cm_8Be : " << e_cm_8Be-m_8Be << " MeV" << G4endl;  
  G4cout << "total lab.px : " << lv_lab_total.px() << G4endl;  
  G4cout << "total lab.py : " << lv_lab_total.py() << G4endl;  
  G4cout << "total lab.pz : " << lv_lab_total.pz() << G4endl;  
  G4cout << "total lab.kinetic : " << lv_lab_total.e()-m_projectile-m_target << G4endl;  
  G4cout << "theta_cm_alpha1 " << theta_cm_alpha1 << G4endl;  
  G4cout << "phi_cm_alpha1 " << phi_cm_alpha1 << G4endl;  
  G4cout << "lv_cm_alpha1.px :        " << lv_cm_alpha1.px() << G4endl;  
  G4cout << "lv_cm_alpha1.py :        " << lv_cm_alpha1.py() << G4endl;  
  G4cout << "lv_cm_alpha1.pz :        " << lv_cm_alpha1.pz() << G4endl;  
  G4cout << "lv_cm_alpha1.kinematic : " << lv_cm_alpha1.e() - m_4He << G4endl;
  G4cout << "lv_cm_8Be.px :        " << lv_cm_8Be.px() << G4endl;  
  G4cout << "lv_cm_8Be.py :        " << lv_cm_8Be.py() << G4endl;  
  G4cout << "lv_cm_8Be.pz :        " << lv_cm_8Be.pz() << G4endl;  
  G4cout << "lv_cm_8Be.kinematic : " << lv_cm_8Be.e() - m_8Be << G4endl; 
  G4cout << "lv_lab_alpha1.px :        " << lv_lab_alpha1.px() << G4endl;  
  G4cout << "lv_lab_alpha1.py :        " << lv_lab_alpha1.py() << G4endl;  
  G4cout << "lv_lab_alpha1.pz :        " << lv_lab_alpha1.pz() << G4endl;  
  G4cout << "lv_lab_alpha1.kinematic : " << lv_lab_alpha1.e() - m_4He << G4endl; 
  G4cout << "lv_lab_alpha1.theta :     " << lv_lab_alpha1.theta()<<G4endl;  
  G4cout << "lv_lab_alpha1.phi :       " << lv_lab_alpha1.phi()<<G4endl;
  G4cout << "lv_lab_8Be.px :        " << lv_lab_8Be.px() << G4endl;  
  G4cout << "lv_lab_8Be.py :        " << lv_lab_8Be.py() << G4endl;  
  G4cout << "lv_lab_8Be.pz :        " << lv_lab_8Be.pz() << G4endl;  
  G4cout << "lv_lab_8Be.kinematic : " << lv_lab_8Be.e() - m_8Be << G4endl; 
  G4cout << "lv_lab_8Be.theta :     " << lv_lab_8Be.theta()<<G4endl;  
  G4cout << "lv_lab_8Be.phi :       " << lv_lab_8Be.phi()<<G4endl;
#endif

  // 
  // 8Be -> alpha+alpha
  G4double s_new = m_8Be*m_8Be;
  
  G4double e_cm_alpha2 = sqrt(s_new)/2.;
  G4double e_cm_alpha3 = sqrt(s_new)/2.;
  G4double p_cm_alpha2 = sqrt(e_cm_alpha2*e_cm_alpha2 - m_4He*m_4He);
  G4double p_cm_alpha3 = sqrt(e_cm_alpha3*e_cm_alpha3 - m_4He*m_4He);
  
  G4ThreeVector dir_cm_alpha2;

  if(flag_seq1 && resonance_type == Resonance165 && H11B165UseExitAngularCorrelation){
    // 165-keV alpha1 branch:
    //   12C*(16.11, 2+) -> alpha + 8Be(2+), L = 2 dominated,
    //   8Be(2+) -> alpha + alpha, L = 2.
    // The sampled direction is in the 8Be rest frame and is correlated with
    // the primary alpha direction in the 12C CM frame.
    dir_cm_alpha2 = H11BAngularDistribution::Sample165Alpha1DWaveDirection(lv_cm_alpha1.vect());
  }
  else if(flag_seq1 && resonance_type == Resonance675 && H11B675UseExitAngularCorrelation){
    // 675-keV alpha1 branch:
    //   12C*(16.57/16.62, 2-) -> alpha + 8Be(2+), coherent L=1/L=3 mixture,
    //   8Be(2+) -> alpha + alpha, L = 2.
    dir_cm_alpha2 = H11BAngularDistribution::Sample675Alpha1L13Direction(lv_cm_alpha1.vect());
  }
  else{
    // alpha0 branch through 8Be(g.s.) or fallback: keep the old isotropic decay.
    dir_cm_alpha2 = H11BAngularDistribution::SampleIsotropicDirection();
  }

  G4double theta_cm_alpha2 = dir_cm_alpha2.theta();
  G4double phi_cm_alpha2 = dir_cm_alpha2.phi();

  G4LorentzVector lv_cm_total_new(0., 0., 0., sqrt(s_new));
  G4LorentzVector lv_cm_alpha2(p_cm_alpha2*dir_cm_alpha2.x(),
                               p_cm_alpha2*dir_cm_alpha2.y(),
                               p_cm_alpha2*dir_cm_alpha2.z(),
                               e_cm_alpha2);
  G4LorentzVector lv_cm_alpha3 = lv_cm_total_new - lv_cm_alpha2;

  G4LorentzVector lv_lab_alpha2 = lv_cm_alpha2;
  lv_lab_alpha2.boost(lv_lab_8Be.boostVector());
  G4LorentzVector lv_lab_alpha3 = lv_cm_alpha3;
  lv_lab_alpha3.boost(lv_lab_8Be.boostVector());

#ifdef ReactionCout
  G4cout << "e_cm_alpha2 : " << e_cm_alpha2 - m_4He << G4endl;  
  G4cout << "e_cm_alpha3 : " << e_cm_alpha3 - m_4He << G4endl;  
  G4cout << "p_cm_alpha2 : " << p_cm_alpha2 << G4endl;  
  G4cout << "p_cm_alpha3 : " << p_cm_alpha3 << G4endl;
  G4cout << "theta_cm_alpha2 " << theta_cm_alpha2 << G4endl;  
  G4cout << "phi_cm_alpha2 " << phi_cm_alpha2 << G4endl;
  G4cout << "lv_cm_alpha2.px : "        << lv_cm_alpha2.px() << G4endl;  
  G4cout << "lv_cm_alpha2.py : "        << lv_cm_alpha2.py() << G4endl;  
  G4cout << "lv_cm_alpha2.pz : "        << lv_cm_alpha2.pz() << G4endl;  
  G4cout << "lv_cm_alpha2.kinematic : " << lv_cm_alpha2.e() - m_4He << G4endl; 
  G4cout << "lv_cm_alpha3.px : "        << lv_cm_alpha3.px() << G4endl;  
  G4cout << "lv_cm_alpha3.py : "        << lv_cm_alpha3.py() << G4endl;  
  G4cout << "lv_cm_alpha3.pz : "        << lv_cm_alpha3.pz() << G4endl;  
  G4cout << "lv_cm_alpha3.kinematic : " << lv_cm_alpha3.e() - m_4He << G4endl; 
  G4cout << "lv_lab_alpha2.px : "        << lv_lab_alpha2.px() << G4endl;  
  G4cout << "lv_lab_alpha2.py : "        << lv_lab_alpha2.py() << G4endl;  
  G4cout << "lv_lab_alpha2.pz : "        << lv_lab_alpha2.pz() << G4endl;  
  G4cout << "lv_lab_alpha2.kinematic : " << lv_lab_alpha2.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha2.theta :     " << lv_lab_alpha2.theta()<<G4endl;  
  G4cout << "lv_lab_alpha2.phi :       " << lv_lab_alpha2.phi()<<G4endl;
  G4cout << "lv_lab_alpha3.px : "        << lv_lab_alpha3.px() << G4endl;  
  G4cout << "lv_lab_alpha3.py : "        << lv_lab_alpha3.py() << G4endl;  
  G4cout << "lv_lab_alpha3.pz : "        << lv_lab_alpha3.pz() << G4endl;  
  G4cout << "lv_lab_alpha3.kinematic : " << lv_lab_alpha3.e() - m_4He << G4endl;
  G4cout << "lv_lab_alpha3.theta :     " << lv_lab_alpha3.theta()<<G4endl;  
  G4cout << "lv_lab_alpha3.phi :       " << lv_lab_alpha3.phi()<<G4endl; 
#endif
  
  G4DynamicParticle *dynamic_product1 = new G4DynamicParticle(product1, lv_lab_alpha1); 
  G4DynamicParticle *dynamic_product2 = new G4DynamicParticle(product2, lv_lab_alpha2); 
  G4DynamicParticle *dynamic_product3 = new G4DynamicParticle(product3, lv_lab_alpha3); 

  // Add the secondaries to the particle change stack
  // Changes being included in process
  theParticleChange.AddSecondary(dynamic_product1);
  theParticleChange.AddSecondary(dynamic_product2);
  theParticleChange.AddSecondary(dynamic_product3);

  // Build the three-alpha center-of-mass frame from the final-state alpha particles.
  // This is the proper frame for the Dalitz plot.
  G4LorentzVector lv_lab_3alpha =
  lv_lab_alpha1 + lv_lab_alpha2 + lv_lab_alpha3;

  G4ThreeVector beta_3alpha = lv_lab_3alpha.boostVector();

  G4LorentzVector lv_3alpha_cm_alpha1 = lv_lab_alpha1;
  G4LorentzVector lv_3alpha_cm_alpha2 = lv_lab_alpha2;
  G4LorentzVector lv_3alpha_cm_alpha3 = lv_lab_alpha3;

  lv_3alpha_cm_alpha1.boost(-beta_3alpha);
  lv_3alpha_cm_alpha2.boost(-beta_3alpha);
  lv_3alpha_cm_alpha3.boost(-beta_3alpha);

  // rootfile
  reaction_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  reaction_data.e_alpha1 = lv_lab_alpha1.e() - m_4He;
  reaction_data.e_alpha2 = lv_lab_alpha2.e() - m_4He;
  reaction_data.e_alpha3 = lv_lab_alpha3.e() - m_4He;
  // three-alpha CM-frame kinetic energies
  reaction_data.e_3alpha_cm_alpha1 = lv_3alpha_cm_alpha1.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha2 = lv_3alpha_cm_alpha2.e() - m_4He;
  reaction_data.e_3alpha_cm_alpha3 = lv_3alpha_cm_alpha3.e() - m_4He;
  reaction_data.theta_lab_alpha1 = lv_lab_alpha1.theta();
  reaction_data.theta_lab_alpha2 = lv_lab_alpha2.theta();
  reaction_data.theta_lab_alpha3 = lv_lab_alpha3.theta();
  reaction_data.phi_lab_alpha1 = lv_lab_alpha1.phi();
  reaction_data.phi_lab_alpha2 = lv_lab_alpha2.phi();
  reaction_data.phi_lab_alpha3 = lv_lab_alpha3.phi();
  reaction_data.ex_8Be = ex_8Be;
  reaction_data.e_8Be = lv_lab_8Be.e() - m_8Be;
  reaction_data.theta_lab_8Be = lv_lab_8Be.theta();
  reaction_data.phi_lab_8Be = lv_lab_8Be.phi();

  reaction_data.x = 0.;
  reaction_data.y = 0.;
  reaction_data.z = 0.;

  if(flag_seq1){
    strcpy(reaction_data.reaction, Form("%s_sequential1", GetResonanceName()));    
  }else if(flag_seq2){
    strcpy(reaction_data.reaction, Form("%s_sequential2", GetResonanceName()));    
  }else if(flag_simu){
    strcpy(reaction_data.reaction, "simultaneous");    
  }else{
    // gamma
  }

  RunAction* run_action = static_cast<RunAction*>(const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
  run_action->GetRootIO()->FillReactionTree(reaction_data);
}


//
G4double H11BReaction::Sample8Be2PlusExcitationEnergy(G4double ex_max)
{
  if(ex_max<=0.) return 0.;

  // The histogram is cached for each kinematic upper limit. The bin width is
  // kept close to 10 keV, similar to the original implementation.
  static std::map<G4int, TH1D*> m_8Be2Plus_hist;
  G4int key = (G4int)(ex_max/keV + 0.5);

  if(m_8Be2Plus_hist.find(key) == m_8Be2Plus_hist.end()){
    G4double bin_width = 10. *keV;
    G4int n_bins = std::max(10, (G4int)(ex_max/bin_width));

    TH1D *h = new TH1D(Form("h_8Be2Plus_exmax_%05dkeV", key), "", n_bins, 0., ex_max);
    h->SetDirectory(nullptr);

    for(G4int i=1;i<=n_bins;i++){
      G4double ex = h->GetBinCenter(i);
      h->SetBinContent(i, Weight8Be2Plus(ex));
    }

    m_8Be2Plus_hist[key] = h;
  }

  return m_8Be2Plus_hist[key]->GetRandom();
}

//
G4double H11BReaction::Weight8Be2Plus(G4double ex_8Be) const
{
  // ex_8Be is relative to the 8Be ground state, which is the quantity added
  // to the 8Be mass in ReactionKinematic(). For the alpha-alpha line shape,
  // the relevant energy is instead the relative energy above the 2-alpha
  // threshold:
  //     Eaa = ex_8Be + Ex8BeGroundAbove2Alpha.
  G4double eaa = ex_8Be + Ex8BeGroundAbove2Alpha;
  if(eaa<=0.) return 0.;

  // Energy-dependent width for 8Be(2+) -> alpha + alpha.
  // Default: full Coulomb penetrability ratio
  //     Gamma(E) = Gamma_R * P_2(E,a) / P_2(E_R,a),
  // where P_L = rho/[F_L(eta,rho)^2 + G_L(eta,rho)^2].
  // If the full-Coulomb flag is disabled, the code falls back to the older
  // threshold approximation Gamma(E)=Gamma_R*(E/E_R)^(5/2).
  G4double gamma_e = Gamma8Be2Plus*std::pow(eaa/Eaa8Be2Plus, 2.5);
  if(H11BUseFullCoulombPenetrabilityFor8Be2Plus){
    const G4double ratio = H11BCoulomb::AlphaAlphaL2Ratio(eaa/keV,
                                                          Eaa8Be2Plus/keV);
    gamma_e = Gamma8Be2Plus*ratio;
  }

  G4double delta_e = eaa - Eaa8Be2Plus;

  return gamma_e/(delta_e*delta_e + gamma_e*gamma_e/4.);
}

//
G4bool H11BReaction::ChooseAlpha0Channel() const
{
  // In this code convention:
  //   seq2 -> 8Be ground state -> alpha0
  //   seq1 -> 8Be excited state -> alpha1
  // The relative alpha0/alpha1 yield is taken directly from the partial
  // widths of the selected resonance.
  G4double probability_alpha0 = GetAlpha0BranchingRatio();
  G4double rand = G4UniformRand();

  if(rand<probability_alpha0) return true;  // alpha0 -> 8Be ground state, seq2
  return false;                            // alpha1 -> 8Be excited state, seq1
}

//
G4double H11BReaction::GetAlpha0BranchingRatio() const
{
  G4double alpha0_width = H11B165Alpha0Width/keV;
  G4double alpha1_width = H11B165Alpha1Width/keV;

  if(resonance_type == Resonance675){
    alpha0_width = H11B675Alpha0Width/keV;
    alpha1_width = H11B675Alpha1Width/keV;
  }

  G4double total_alpha_width = alpha0_width + alpha1_width;
  if(total_alpha_width<=0.) return 0.;

  return alpha0_width/total_alpha_width;
}

//
const char* H11BReaction::GetResonanceName() const
{
  if(resonance_type == Resonance675) return "p11B_675";
  return "p11B_165";
}
