#include "JUNAReaction.hh"

#include <TGenPhaseSpace.h>
#include <TMath.h>
//#include <TTree.h>
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4ParticleMomentum.hh"
#include "G4IonTable.hh"
#include "G4IonTable.hh"
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
//#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4AnalysisManager.hh"
#include "TF1.h"
#include "TGraph.h"
//#include "func.C" 
Double_t ff(Double_t *x, Double_t *) { return (1 + 2/7*(3*x[0]*x[0]-1)/2 - 9/7*(35*x[0]*x[0]*x[0]*x[0] - 30*x[0]*x[0] + 3)/8); }
Double_t axisChange(Double_t vx,Double_t vy,Double_t vz,Double_t v0,Double_t* v,Double_t theta,Double_t phi,Double_t* thetal,Double_t* phil)
{
Double_t vxl = vz*sin(theta)*cos(phi) + vx*cos(theta)*cos(phi) + vy*cos(phi+pi/2) + v0*sin(theta)*cos(phi);
Double_t vyl = vz*sin(theta)*sin(phi) + vx*cos(theta)*sin(phi) + vy*sin(phi+pi/2) + v0*sin(theta)*sin(phi);
Double_t vzl = vz*cos(theta) - vx*sin(theta) + v0*cos(theta);
*v = sqrt(vxl*vxl + vyl*vyl + vzl*vzl);
*thetal = acos(vzl/ *v);

if(vxl>0 && vyl>0) *phil = atan(vyl/vxl);
if(vxl<0 && vyl>0) *phil = atan(vyl/vxl) + pi;
if(vxl<0 && vyl<0) *phil = atan(vyl/vxl) + pi;
if(vxl>0 && vyl<0) *phil = atan(vyl/vxl) + 2*pi;
return 0;
}

Double_t GetP(Double_t E1)
{
   Double_t p;
   Double_t a[27],b[27];
   a[0] = 51.3; b[0] = 2.9/213;
   a[1] = 59.6; b[1] = 2.57/209;
   a[2] = 67.8; b[2] = 2.12/205;
   a[3] = 72.6; b[3] = 2.8/201;
   a[4] = 77.6; b[4] = 2.12/218;
   a[5] = 84.3; b[5] = 2.34/213;
   a[6] = 93.1; b[6] = 2.4/235;
   a[7] = 102.1;b[7] = 2.3/234;
   a[8] = 111.2;b[8] = 2.62/266;
   a[9] = 120.3;b[9] = 2.82/255;
   a[10] = 129.3;b[10] = 3.67/290;
   a[11] = 135.9;b[11] = 5/328;
   a[12] = 162.9;b[12] = 5.2/392;
   a[13] = 171.9;b[13] = 3.2/304;
   a[14] = 181; b[14] = 2.54/284;
   a[15] = 190;b[15] = 2.16/262;
   a[16] = 199;b[16] = 2.08/265;
   a[17] = 208;b[17] = 1.87/241;
   a[18] = 221;b[18] = 1.91/266;
   a[19] = 235;b[19] = 2.1/270;
   a[20] = 244;b[20] = 1.98/269;
   a[21] = 253;b[21] = 1.98/273;
   a[22] = 266;b[22] = 1.9/271;
   a[23] = 284;b[23] = 1.77/282;
   a[24] = 302;b[24] = 1.79/280;
   a[25] = 316;b[25] = 1.79/290;
   a[26] = 325;b[26] = 1.77/315;
   int i = 0;
   if (E1<a[0]) p = b[0];
     else if (E1>a[26]) p = b[26];
       else {
         for (i=1;i<27;i++)
         {
            if (E1<a[i]) 
            {
                p = (b[i]+b[i-1])/2; 
                break;
            }
         }
       }      
   return p/(p+1);
}


JUNAReaction::JUNAReaction(): G4HadronicInteraction()
{

  SetMinEnergy(0. * CLHEP::keV);
  SetMaxEnergy(100. * CLHEP::MeV);
  isBlocked = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
JUNAReaction::~JUNAReaction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadFinalState *JUNAReaction::ApplyYourself(const G4HadProjectile &projectile, G4Nucleus &targetNucleus)
{
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill);

  G4LorentzVector proLonVector = projectile.Get4Momentum();
  G4int targetZ = targetNucleus.GetZ_asInt();
  const G4int targetA = targetNucleus.GetA_asInt();
  const G4int projectA = projectile.GetDefinition()->GetBaryonNumber();

  if(projectA==1 && targetA==11){
    G4ParticleDefinition *target = G4IonTable::GetIonTable()->GetIon(targetZ, targetA, 0. * CLHEP::keV);
  //  G4ParticleDefinition *product0 = G4Gamma::GammaDefinition();
 //   G4ParticleDefinition *product1 = G4IonTable::GetIonTable()->GetIon(targetZ+1, targetA+1, 0. * keV);
    G4ParticleDefinition *product0 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition *product1 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    G4ParticleDefinition *product2 = G4IonTable::GetIonTable()->GetIon(2, 4, 0. * CLHEP::keV);
    ReactionKinematic(projectile, target, product0, product1, product2);
  }

  return &theParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void JUNAReaction::ReactionKinematic(const G4HadProjectile &projectile, G4ParticleDefinition *targetNucleus, G4ParticleDefinition *product0, G4ParticleDefinition *product1, G4ParticleDefinition *product2)
{
  G4double Ek0 = projectile.GetKineticEnergy(); //MeV
  G4double Ek1 = Ek0/12;
  
  G4double EBe = 3.03 + (G4UniformRand()-0.5)*1.513;
  G4double Pos0 = GetP(Ek0/12*11);
  if (G4UniformRand() < Pos0)  EBe = 0;
  
  G4double Ex1 = (Ek0/12*11 + (15.9567-7.3666-EBe)) ;
  G4double Ex3 = 0.0919+EBe;
  
  
  G4ThreeVector DC = projectile.GetMomentumDirection();
  G4double theta0 = DC.theta();
  G4double phi0 = DC.phi();
  if(phi0<0) phi0 = phi0 +2*pi;
  //G4cout<<" Ek0 = "<<Ek0<<" Theta: "<<Tproj<<" Phi: "<<Pproj<<G4endl;
  
  
  G4double mProject = 1*1.66054*1e-27; // kg
  G4double mTarget = 11*1.66054*1e-27;
  G4double mC = 12*1.66054*1e-27;
  G4double mAlpha = 4*1.66054*1e-27;
  G4double mBe = 8*1.66054*1e-27;
  
  G4double MeV2J = 1.60217646*1e-13;
  
  G4double theta1 = (G4UniformRand()-0.5)*2;
  theta1 = acos(theta1);
  G4double theta2 = pi - theta1;
  
  G4double VC = sqrt(2*Ek1*MeV2J/mC);
  Double_t V1 = sqrt(Ex1*MeV2J/(0.5*mAlpha + mAlpha*mAlpha/2/mBe));
  Double_t V2 = mAlpha/mBe*V1;
  G4double phi1 = G4UniformRand()*2*pi;
  G4double phi2 = phi1 + pi;

  if(phi2 > 2*pi) phi2 = phi2 - 2*pi;
  
  G4double v1x = V1*sin(theta1)*cos(phi1);
  G4double v1y = V1*sin(theta1)*sin(phi1);
  G4double v1z = V1*cos(theta1);
  G4double v2x = V2*sin(theta2)*cos(phi2);
  G4double v2y = V2*sin(theta2)*sin(phi2);
  G4double v2z = V2*cos(theta2);
  Double_t theta1l,phi1l,theta2l,phi2l;
  axisChange(v1x,v1y,v1z,VC,&V1,theta0,phi0,&theta1l,&phi1l);
  axisChange(v2x,v2y,v2z,VC,&V2,theta0,phi0,&theta2l,&phi2l);
  
  TF1 *fu = new TF1("fu", ff, -1, 1);
  //G4double theta3 = acos(fu->GetRandom());
  G4double theta3 = acos((G4UniformRand()-0.5)*2);
  G4double theta4 = pi - theta3;
  G4double V3 = sqrt(Ex3*MeV2J/(0.5*mAlpha + mAlpha*mAlpha/2/mAlpha));
  G4double V4 = V3;
  
  G4double phi3 = G4UniformRand()*2*pi;
  G4double phi4 = phi3 + pi;
    
  if(phi2 > 2*pi) phi2 = phi2 - 2*pi;
  G4double v3x = V3*sin(theta3)*cos(phi3);
  G4double v3y = V3*sin(theta3)*sin(phi3);
  G4double v3z = V3*cos(theta3);
  G4double v4x = V4*sin(theta4)*cos(phi4);
  G4double v4y = V4*sin(theta4)*sin(phi4);
  G4double v4z = V4*cos(theta4);
  Double_t theta3l,phi3l,theta4l,phi4l;
  axisChange(v3x,v3y,v3z,V2,&V3,theta2l,phi2l,&theta3l,&phi3l);
  axisChange(v4x,v4y,v4z,V2,&V4,theta2l,phi2l,&theta4l,&phi4l);

  G4double Pconsx = mC*VC*sin(theta0)*cos(phi0) - mAlpha*V1*sin(theta1l)*cos(phi1l) - mAlpha*V3*sin(theta3l)*cos(phi3l) - mAlpha*V4*sin(theta4l)*cos(phi4l);
  G4double Pconsy = mC*VC*sin(theta0)*sin(phi0) - mAlpha*V1*sin(theta1l)*sin(phi1l) - mAlpha*V3*sin(theta3l)*sin(phi3l) - mAlpha*V4*sin(theta4l)*sin(phi4l);
  G4double Pconsz = mC*VC*cos(theta0) - mAlpha*V1*cos(theta1l) - mAlpha*V3*cos(theta3l) - mAlpha*V4*cos(theta4l);
  G4double Econs = (mAlpha*V1*V1/2 + mAlpha*V3*V3/2 + mAlpha*V4*V4/2)  -Ek1  - Ex1 - Ex3;
  
  G4double alphaMass = product0->GetPDGMass()*CLHEP::MeV;
  G4double p1x = mAlpha*V1*sin(theta1l)*cos(phi1l)*299792458/(1.60271646*1e-13); //to MeV/c
  G4double p1y = mAlpha*V1*sin(theta1l)*sin(phi1l)*299792458/(1.60271646*1e-13);
  G4double p1z = mAlpha*V1*cos(theta1l)*299792458/(1.60271646*1e-13);
  G4double e1 =  mAlpha*V1*V1/2/(1.60271646*1e-13) + alphaMass;
  G4double p3x = mAlpha*V3*sin(theta3l)*cos(phi3l)*299792458/(1.60271646*1e-13); //to MeV/c
  G4double p3y = mAlpha*V3*sin(theta3l)*sin(phi3l)*299792458/(1.60271646*1e-13);
  G4double p3z = mAlpha*V3*cos(theta3l)*299792458/(1.60271646*1e-13);
  G4double e3 =  mAlpha*V3*V3/2/(1.60271646*1e-13) + alphaMass;
  G4double p4x = mAlpha*V4*sin(theta4l)*cos(phi4l)*299792458/(1.60271646*1e-13); //to MeV/c
  G4double p4y = mAlpha*V4*sin(theta4l)*sin(phi4l)*299792458/(1.60271646*1e-13);
  G4double p4z = mAlpha*V4*cos(theta4l)*299792458/(1.60271646*1e-13);
  G4double e4 =  mAlpha*V4*V4/2/(1.60271646*1e-13) + alphaMass;
  /*
  G4cout<<"p1x "<<p1x<<"p1y "<<p1y<<"p1z "<<p1z<<"e1 "<<e1<<G4endl;
  G4cout<<"V1 "<<V1<<"V2 "<<V2<<"V3 "<<V3<<"V4 "<<V4<<G4endl;
  G4cout<<" Px "<<Pconsx<<" Py "<<Pconsy<< " Pz "<<Pconsz << " EE "<<Econs<< G4endl;
  */
  
G4LorentzVector product0LonVector, product1LonVector, product2LonVector;

  product0LonVector.setPx(p1x);
  product0LonVector.setPy(p1y);
  product0LonVector.setPz(p1z);
  product0LonVector.setE(e1);
  product1LonVector.setPx(p3x);
  product1LonVector.setPy(p3y);
  product1LonVector.setPz(p3z);
  product1LonVector.setE(e3);
  product2LonVector.setPx(p4x);
  product2LonVector.setPy(p4y);
  product2LonVector.setPz(p4z);
  product2LonVector.setE(e4);

  G4DynamicParticle *product0Dynamic, *product1Dynamic,*product2Dynamic;
  product0Dynamic = new G4DynamicParticle(product0, product0LonVector);
  product1Dynamic = new G4DynamicParticle(product1, product1LonVector);
  product2Dynamic = new G4DynamicParticle(product2, product2LonVector);

  // Add the secondaries to the particle change stack
  // Changes being included in process
  theParticleChange.AddSecondary(product0Dynamic);
  theParticleChange.AddSecondary(product1Dynamic);
  theParticleChange.AddSecondary(product2Dynamic);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager -> FillNtupleDColumn(1,0,e1-alphaMass);
  analysisManager -> FillNtupleDColumn(1,1,theta1l);
  analysisManager -> FillNtupleDColumn(1,2,phi1l);
  analysisManager -> FillNtupleDColumn(1,3,e3-alphaMass);
  analysisManager -> FillNtupleDColumn(1,4,theta3l);
  analysisManager -> FillNtupleDColumn(1,5,phi3l);
  analysisManager -> FillNtupleDColumn(1,6,e4-alphaMass);
  analysisManager -> FillNtupleDColumn(1,7,theta4l);
  analysisManager -> FillNtupleDColumn(1,8,phi4l);
  analysisManager -> AddNtupleRow(1);
}
