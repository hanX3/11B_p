#ifndef H11BConfig_H
#define H11BConfig_H 1

#include "Constants.hh"
#include "G4String.hh"
#include "globals.hh"

#include <utility>

namespace H11BConfig {
G4bool GetEnable165PrimaryAngularDistribution();
void SetEnable165PrimaryAngularDistribution(G4bool enabled);
void Set165PrimaryAngularA1(G4double value);
void Set165PrimaryAngularA2(G4double value);
G4double Get165PrimaryAngularA1();
G4double Get165PrimaryAngularA2();
std::pair<G4double, G4double> Get165PrimaryA1A2(G4bool alpha1_branch);

G4bool GetEnable675PrimaryAngularDistribution();
void SetEnable675PrimaryAngularDistribution(G4bool enabled);
void Set675PrimaryAngularA1(G4double value);
void Set675PrimaryAngularA2(G4double value);
G4double Get675PrimaryAngularA1();
G4double Get675PrimaryAngularA2();

G4bool GetEnable165Alpha1SecondaryAngularCorrelation();
void SetEnable165Alpha1SecondaryAngularCorrelation(G4bool enabled);
const char* Get165Alpha1SecondaryAngularCorrelationName();
void Set165Alpha1SecondaryA2(G4double value);
void Set165Alpha1SecondaryA4(G4double value);
G4double Get165Alpha1SecondaryA2();
G4double Get165Alpha1SecondaryA4();

G4bool GetEnable675Alpha1SecondaryAngularCorrelation();
void SetEnable675Alpha1SecondaryAngularCorrelation(G4bool enabled);
const char* Get675Alpha1SecondaryAngularCorrelationName();
void Set675Alpha1SecondaryA2(G4double value);
void Set675Alpha1SecondaryA4(G4double value);
G4double Get675Alpha1SecondaryA2();
G4double Get675Alpha1SecondaryA4();

G4bool GetEnable165Gamma0AngularDistribution();
void SetEnable165Gamma0AngularDistribution(G4bool enabled);
const char* Get165Gamma0AngularDistributionName();
H11BGammaAngularMode Get165Gamma0AngularModeForOutput();
void Set165Gamma0AngularA1(G4double value);
void Set165Gamma0AngularA2(G4double value);
G4double Get165Gamma0AngularA1();
G4double Get165Gamma0AngularA2();

G4bool GetEnable675GammaAngularDistribution();
void SetEnable675GammaAngularDistribution(G4bool enabled);
const char* Get675GammaAngularDistributionName();
H11BGammaAngularMode Get675GammaAngularModeForOutput();
void Set675GammaAngularA1(G4double value);
void Set675GammaAngularA2(G4double value);
G4double Get675GammaAngularA1();
G4double Get675GammaAngularA2();

G4bool Get165GammaCaptureEnabled();
void Set165GammaCaptureEnabled(G4bool enabled);

G4bool Get675GammaCaptureEnabled();
void Set675GammaCaptureEnabled(G4bool enabled);

G4double GetGammaBiasFactor();
void SetGammaBiasFactor(G4double factor);

} // namespace H11BConfig

#endif
