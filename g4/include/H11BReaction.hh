#ifndef H11BReaction_H
#define H11BReaction_H 1

#include "RunAction.hh"
#include "G4String.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"

#include "DataStructure.hh"
#include "H11BCrossSection.hh"
#include "RootIO.hh"
#include "Constants.hh"

#include <memory>

class G4GenericMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class H11BReaction : public G4HadronicInteraction
{
public:
  enum class ReactionChannel
  {
    Resonance165 = 0,
    Resonance675 = 1,
    Background3Alpha = 2,
    GammaCapture12C = 3
  };

  enum class BackgroundMode
  {
    PhaseSpace = 2
  };

  enum ResonanceType
  {
    Resonance165,
    Resonance675
  };

public:
  H11BReaction();
  ~H11BReaction() override;

public:
  G4HadFinalState* ApplyYourself(const G4HadProjectile& projectile, G4Nucleus& target) override;
  void ReactionKinematic(const G4HadProjectile& projectile, G4ParticleDefinition* target, G4ParticleDefinition* product1, G4ParticleDefinition* product2, G4ParticleDefinition* product3);
  static BackgroundMode GetBackgroundMode();
  static void SetBackgroundMode(BackgroundMode mode);
  static void SetBackgroundMode(const G4String& mode);
  static const char* GetBackgroundModeName();
  void SetEvaluatedCrossSectionModeCommand(const G4String& mode);
  void SetCrossSectionBiasFactorCommand(G4double factor);
  void SetBackgroundBiasFactorCommand(G4double factor);
  void SetBackgroundModeCommand(const G4String& mode);
  void SetEnable675PrimaryAngularDistributionCommand(G4bool enabled);
  void Set675PrimaryAngularA1Command(G4double value);
  void Set675PrimaryAngularA2Command(G4double value);
  void SetEnable675Alpha1SecondaryAngularCorrelationCommand(G4bool enabled);
  void Set675Alpha1SecondaryA2Command(G4double value);
  void Set675Alpha1SecondaryA4Command(G4double value);
  void SetEnable165PrimaryAngularDistributionCommand(G4bool enabled);
  void Set165PrimaryAngularA1Command(G4double value);
  void Set165PrimaryAngularA2Command(G4double value);
  void SetEnable165Alpha1SecondaryAngularCorrelationCommand(G4bool enabled);
  void Set165Alpha1SecondaryA2Command(G4double value);
  void Set165Alpha1SecondaryA4Command(G4double value);
  void SetEnable165Gamma0AngularDistributionCommand(G4bool enabled);
  void Set165Gamma0AngularA1Command(G4double value);
  void Set165Gamma0AngularA2Command(G4double value);
  void SetEnable675GammaAngularDistributionCommand(G4bool enabled);
  void Set675GammaAngularA1Command(G4double value);
  void Set675GammaAngularA2Command(G4double value);
  void Set165GammaCaptureEnabledCommand(G4bool enabled);
  void Set675GammaCaptureEnabledCommand(G4bool enabled);
  void SetGammaBiasFactorCommand(G4double factor);
  void PrintConfigCommand();

private:
  H11BReactionData reaction_data;
  ResonanceType resonance_type = Resonance165;
  ReactionChannel reaction_channel = ReactionChannel::Resonance165;
  H11BGammaResonance gamma_resonance = H11BGammaResonance::None;
  H11BGammaBranch gamma_branch = H11BGammaBranch::None;
  H11BCrossSectionComponents selected_components;
  std::unique_ptr<G4GenericMessenger> messenger;

private:
  G4bool SelectReactionChannel(G4double kinetic_energy_lab);
  void DefineCommands();
  void GenerateThreeBodyPhaseSpace(const G4HadProjectile& projectile, G4ParticleDefinition* target, G4ParticleDefinition* product1, G4ParticleDefinition* product2, G4ParticleDefinition* product3);
  void GenerateGammaCapture(const G4HadProjectile& projectile, G4ParticleDefinition* target);
  void Generate165Gamma0(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle, const G4ThreeVector& beam_axis_cm);
  void Generate165Gamma1Cascade(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle);
  void Generate675GammaLine(const G4LorentzVector& initial_cm, const G4LorentzVector& initial_lab, G4ParticleDefinition* c12_ground, G4ParticleDefinition* gamma_particle, const G4ThreeVector& beam_axis_cm);
  G4double SamplePhaseSpaceM23(G4double sqrt_s, G4double m_alpha) const;
  G4double Sample8Be2PlusExcitationEnergy(G4double ex_max);
  G4double Weight8Be2Plus(G4double ex_8Be, G4double ex_max) const;
  G4double WeightAlpha8BeFirstBreakup(G4double e_alpha8Be, G4double e_alpha8Be_ref) const;
  G4bool ChooseAlpha0Channel() const;
  G4double GetAlpha0BranchingRatio() const;
  const char* GetResonanceName() const;
};

#endif
