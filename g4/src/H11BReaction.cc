#include "H11BReaction.hh"

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

#include "TF1.h"
#include "TGraph.h"
#include <TGenPhaseSpace.h>
#include <TMath.h>

#include "Constants.hh"

#include "Randomize.hh"

//
TSpline3* H11BReaction::sp3_seq2over1 = nullptr;

//
H11BReaction::H11BReaction()
: G4HadronicInteraction()
{
  if(!sp3_seq2over1) SetSeq2Over1Spline3();

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

  //
  // p+11B -> alpha+8Be
  G4double ex_8Be = 0.;
  G4bool flag_seq1 = 0, flag_seq2 = 0, flag_simu = 0;
  if(!ChooseSeq(kinetic_lab_projectile*m_target/(m_projectile+m_target) / keV)) flag_seq1 = 1;
  else flag_seq2 = 1;
  if(flag_seq1) ex_8Be = GetExBreitW(Ex8Be, Ex8BeGamma);

  m_8Be += ex_8Be;
  
  G4double s = m_projectile*m_projectile + m_target*m_target + 2*m_target*(m_projectile + kinetic_lab_projectile);
  
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
  
  G4double theta_cm_alpha1 = std::acos((G4UniformRand()-0.5)*2.);
  G4double phi_cm_alpha1 = G4UniformRand()*2*TMath::Pi();

  G4LorentzVector lv_cm_alpha1(p_cm_alpha1*std::sin(theta_cm_alpha1)*std::cos(phi_cm_alpha1), 
                               p_cm_alpha1*std::sin(theta_cm_alpha1)*sin(phi_cm_alpha1), 
                               p_cm_alpha1*std::cos(theta_cm_alpha1), e_cm_alpha1);
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
  
  G4double theta_cm_alpha2 = std::acos((G4UniformRand()-0.5)*2.);
  G4double phi_cm_alpha2 = G4UniformRand()*2*TMath::Pi();

  G4LorentzVector lv_cm_total_new(0., 0., 0., sqrt(s_new));
  G4LorentzVector lv_cm_alpha2(p_cm_alpha2*std::sin(theta_cm_alpha2)*std::cos(phi_cm_alpha2), 
                               p_cm_alpha2*std::sin(theta_cm_alpha2)*std::sin(phi_cm_alpha2), 
                               p_cm_alpha2*std::cos(theta_cm_alpha2), e_cm_alpha2);
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

  // rootfile
  reaction_data.event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  reaction_data.e_alpha1 = lv_lab_alpha1.e() - m_4He;
  reaction_data.e_alpha2 = lv_lab_alpha2.e() - m_4He;
  reaction_data.e_alpha3 = lv_lab_alpha3.e() - m_4He;
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
    strcpy(reaction_data.reaction, "sequential1");    
  }else if(flag_seq2){
    strcpy(reaction_data.reaction, "sequential2");    
  }else if(flag_simu){
    strcpy(reaction_data.reaction, "simultaneous");    
  }else{
    // gamma
  }

  RunAction* run_action = static_cast<RunAction*>(const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
  run_action->GetRootIO()->FillReactionTree(reaction_data);
}

//
G4double H11BReaction::GetExBreitW(G4double ex_aver, G4double gamma)
{
  static std::map<std::pair<G4int,G4int>, TH1D*> m_ex_gamma_hist;
  auto key = std::make_pair((G4int)(ex_aver*1000), (G4int)(gamma*1000));
  
  if(m_ex_gamma_hist.find(key) == m_ex_gamma_hist.end()){
     TH1D *h = new TH1D(Form("h_ex_%03dkeV_gamma_%03dkeV", key.first, key.second), "", 2*ex_aver*100, 0, 2*ex_aver);
     
    for(int i=0;i<2*ex_aver*100;i++){
      h->Fill(i/100., gamma*gamma/4./((i/100.-ex_aver)*(i/100.-ex_aver)+gamma*gamma/4.));
    }
    m_ex_gamma_hist[key] = h;
  }

  return m_ex_gamma_hist[key]->GetRandom();
}

//
void H11BReaction::SetSeq2Over1Spline3()
{
  if (sp3_seq2over1) return;

  G4double e_cm[27] = {51.3, 59.6, 67.8, 72.6, 77.6, 84.3, 93.1, 102.1, 111.2, 120.3, 
                       129.3, 135.9, 162.9, 171.9, 181.0, 190.0, 199.0, 208.0, 221.0, 235.0,
                       244.0, 253.0, 266.0, 284.0, 302.0, 316.0, 325.0};

  G4double seq2over1[27] = {2.9/213., 2.57/209., 2.12/205., 2.8/201., 2.12/218.,
                               2.34/213., 2.4/235., 2.3/234., 2.62/266., 2.82/255.,
                               3.67/290., 5.0/328., 5.2/392., 3.2/304., 2.54/284.,
                               2.16/262., 2.08/265., 1.87/241., 1.91/266., 2.1/270.,
                               1.98/269., 1.98/273., 1.9/271., 1.77/282., 1.79/280.,
                               1.79/290., 1.77/315.};

  TGraph *gr_seq2over1 = new TGraph(27, e_cm, seq2over1);
  sp3_seq2over1 = new TSpline3("seq2over1", gr_seq2over1);
}

//
G4bool H11BReaction::ChooseSeq(G4double e_cm)
{
  G4double val = sp3_seq2over1->Eval(e_cm);

  G4double probability = val/(val+1.);

  G4double rand = G4UniformRand();
  // G4cout << "val " << val << G4endl;
  // G4cout << "e_cm " << e_cm << " probability " << probability << " rand " << rand << G4endl;

  if(rand<probability) return 1; // ex_8Be=0, seq2
  else return 0; // ex_8Be>0, seq1
}
