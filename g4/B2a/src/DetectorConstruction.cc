//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "StripSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Element.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

#include "G4StepLimiterPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


DetectorConstruction::DetectorConstruction()
{

  fLogicStrip = new G4LogicalVolume*[30];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete [] fLogicStrip;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Lead defined using NIST Manager
  //fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Pt");

  // Xenon gas defined using NIST Manager
 // G4Element *B = new G4Element("Boron","B",5,10.81*g/mole);
 // G4Material *B_M = new G4Material("B_M",2.37*g/cm3,1);
 // B_M -> AddElement(B,1);
  fStripMaterial = nistManager->FindOrBuildMaterial("G4_Si");
  //-----------------------------name--Z----A----
  G4Isotope* B11 = new G4Isotope("B11",5 , 11, 11*g/mole);
  G4Isotope* B10 = new G4Isotope("B10",5 , 10, 10*g/mole);
  //-----------------------------------name-----symbol--ncomponets
  G4Element* TargetEl = new G4Element("TargetEl", "B",    2);
  TargetEl -> AddIsotope(B11,99*perCent);
  TargetEl -> AddIsotope(B10,1*perCent);
  //------------------------------------name---------density----ncomponets
  fTargetMaterial= new G4Material("TargetB",  1.404*g/cm3,     1);
  fTargetMaterial -> AddElement(TargetEl,1);
  
  
  
    // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* Nothing  = nistManager->FindOrBuildMaterial("G4_Galactic");

  // Sizes of the principal geometrical components (solids)
  G4double worldLength = 40*cm;
  G4double targetThick =  0.239*um; // full length of Target
  G4double targetLong =  20*mm;
  G4double targetShort = 10*mm;
  G4ThreeVector targetPos = G4ThreeVector(0,0,targetThick/2);

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4Box* worldS
    = new G4Box("world",                                    //its name
                worldLength/2,worldLength/2,worldLength/2); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,   //its solid
                 Nothing,      //its material
                 "World"); //its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 0); // checking overlaps

  // Target

  G4Box* targetS
    = new G4Box("target",targetShort/2,targetLong/2,targetThick/2);
   fLogicTarget
    = new G4LogicalVolume(targetS, fTargetMaterial,"Target_LV",0,0,0);
  new G4PVPlacement(0,               // no rotation
                    targetPos,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target_PV",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    0); // checking overlaps
  // Tracker

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  worldLV      ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  
  G4Tubs* BigSi
        = new G4Tubs("BigSi", 0, 1.4*cm, 75*um, 0 ,360);
        
  G4LogicalVolume* BigSiLogic = new G4LogicalVolume(BigSi,fStripMaterial,"BigSi_LV",0,0,0);
  auto rotation = new G4RotationMatrix();
 // rotation->rotateY(45*deg);
 rotation->rotateY(90*deg);
  G4ThreeVector BigSiPos;
  //BigSiPos = G4ThreeVector(1.414*5*cm,0,-1.414*5*cm);
  BigSiPos = G4ThreeVector(4*cm,0,0);
  new G4PVPlacement(rotation,
                      BigSiPos,
                      BigSiLogic,
                      "BigSi_PV",
                      worldLV,
                      false,
                      30,
                      0);

  for (G4int copyNo=0; copyNo<30; copyNo++) {

      G4Box* strip
        = new G4Box("Strip",1*mm,75*um,5*mm);

      fLogicStrip[copyNo] =
              new G4LogicalVolume(strip,fStripMaterial,"Strip_LV",0,0,0);

      fLogicStrip[copyNo]->SetVisAttributes(chamberVisAtt);
      G4ThreeVector stripPos;
      if (copyNo<15)
      stripPos = G4ThreeVector((15.75 - copyNo*2.25)*mm,40*mm,-5*mm);
      if (copyNo>14)
      stripPos = G4ThreeVector((15.75 - (copyNo-15)*2.25)*mm,-40*mm,-5*mm);
      new G4PVPlacement(0,                            // no rotation
                        stripPos, // at (x,y,z)
                        fLogicStrip[copyNo],        // its logical volume
                        "Strip_PV",                 // its name
                        worldLV,                    // its mother  volume
                        false,                        // no boolean operations
                        copyNo,                       // copy number
                        0);              // checking overlaps

  }

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  //G4double maxStep = 0.5*chamberWidth;
  //fStepLimit = new G4UserLimits(maxStep);
  //trackerLV->SetUserLimits(fStepLimit);

  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world
  G4UserLimits* stepLimit = new G4UserLimits();
  stepLimit ->SetMaxAllowedStep(0.02*um);
 
  fLogicTarget -> SetUserLimits(stepLimit);
 
  fScoringVolume = fLogicTarget;
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerStripSDname = "/StripSD";
  StripSD* aStripSD = new StripSD(trackerStripSDname,
                                            "StripHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aStripSD);
  // Setting aTrackerSD to all logical volumes with the same name
  // of "Chamber_LV".
  SetSensitiveDetector("Strip_LV", aStripSD, true);
  SetSensitiveDetector("BigSi_LV", aStripSD, true);



  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
 // G4ThreeVector fieldValue = G4ThreeVector();
 // fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  //fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
 // G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

