#include "H11BCrossSection.hh"

//
H11BCrossSection::H11BCrossSection()
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
  
  const G4double EG = 22.589e3;
  
  
  if(aParA==1 && Z==5 && A==11){
    G4double ecm = GetEcmValue(aParA, A, aPar->GetKineticEnergy() / keV); //unit:keV
    if(ecm>=1 && ecm<=400){
      crossSection = 1e-24*(1.97e5 + 0.24e3*ecm + 2.31e-1*pow(ecm,2) + 1.82e7/(pow((ecm-148),2)+pow(2.35,2))) / ecm*exp(-sqrt(EG/ecm));  
    } //unit cm2
    else if(ecm>400 && ecm<=642){
       crossSection = 1e-24*(3.30e5 + 66.1e3*(ecm/100-4) - 20.3e3*pow((ecm/100-4),2) - 1.58e3*(pow((ecm/100-4),5))) / ecm*exp(-sqrt(EG/ecm));
    }
    else if(ecm>642 && ecm<=3500){
       crossSection = 1e-24*(4.38e3 + 2.57e9/(pow(ecm-581.3,2) + pow(85.7,2)) + 5.67e8/(pow(ecm-1083,2) + pow(234,2)) + 1.34e8/(pow(ecm-2405,2) + pow(138,2)) + 5.68e8/(pow(ecm-3344,2) + pow(309,2))) / ecm*exp(-sqrt(EG/ecm));
    }
    else{
      crossSection = 0.;
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
