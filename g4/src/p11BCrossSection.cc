#include "p11BCrossSection.hh"


p11BCrossSection::p11BCrossSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
p11BCrossSection::~p11BCrossSection()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double p11BCrossSection::GetIsoCrossSection(const G4DynamicParticle *aPar, G4int Z, G4int A, const G4Isotope *, const G4Element *, const G4Material *)
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
    crossSection *= 1e10;

    /*
    if(ecm>=0.2 && ecm<0.25){
      sfactor = 0.1849*ecm*ecm-0.01524*ecm+0.006219;
    }
    else if(ecm>=0.25 && ecm<0.3){
      sfactor = 3.505*pow(ecm,3)-2.065*pow(ecm,2)+0.4148*ecm-0.02436;
    }
    else if(ecm>=0.3 && ecm<0.35){
      sfactor = 35.18*pow(ecm,3)-30.75*pow(ecm,2)+9.076*ecm-0.8964;
    }
    else if(ecm>=0.3 && ecm<0.4){
      sfactor = 1.204*pow(10,-16)*exp(84.78*ecm)+1.115*pow(10,-5)*exp(21.63*ecm);
    }
    else{
      crossSection = 0.;
    }
    
    G4double eta = 0.1575*1.*6*sqrt((1.*12.)/(1.+12.)/ecm);
    crossSection = sfactor*exp(-2.*3.1415926*eta)/ecm * barn;
    crossSection *= 1e10;
    */
  }
  else{
      crossSection = 0.;
  }
  crossSection = crossSection*cm2;
 // G4cout<<crossSection<<G4endl;
  return crossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double p11BCrossSection::GetEcmValue(G4double projectA, G4double targetA, G4double elab)
{
  return targetA/(projectA+targetA)*elab;
}
