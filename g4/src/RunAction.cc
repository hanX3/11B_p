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
/// \file RunAction.cc
/// \brief Implementation of the B2::RunAction class

#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(100000);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager -> SetVerboseLevel(1);
  analysisManager -> CreateNtuple("HB","Energy CopyNb TrackID");
  analysisManager -> CreateNtupleDColumn(0,"Energy");
  analysisManager -> CreateNtupleDColumn(0,"CopyNb");
  analysisManager -> CreateNtupleDColumn(0,"TrackID");
  analysisManager -> FinishNtuple(0);
  analysisManager -> CreateNtuple("AllReaction","E1 theta1 phi1 E2 theta2 phi2 E3 theta3 phi3 ExBe");
  analysisManager -> CreateNtupleDColumn(1,"E1");
  analysisManager -> CreateNtupleDColumn(1,"theta1");
  analysisManager -> CreateNtupleDColumn(1,"phi1");

  analysisManager -> CreateNtupleDColumn(1,"E2");
  analysisManager -> CreateNtupleDColumn(1,"theta2");
  analysisManager -> CreateNtupleDColumn(1,"phi2");

  analysisManager -> CreateNtupleDColumn(1,"E3");
  analysisManager -> CreateNtupleDColumn(1,"theta3");
  analysisManager -> CreateNtupleDColumn(1,"phi3");
  
  analysisManager -> CreateNtupleDColumn(1,"EBe");
  analysisManager -> CreateNtupleDColumn(1,"thetaBe");
  analysisManager -> CreateNtupleDColumn(1,"phiBe");
  
  analysisManager -> CreateNtupleDColumn(1,"ExBe");

  analysisManager -> FinishNtuple(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  auto analysisManager = G4AnalysisManager::Instance();
  G4String fileName = "../data/HBresult.root";
  analysisManager -> OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
auto analysisManager = G4AnalysisManager::Instance();
   analysisManager -> Write();
   analysisManager -> CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


