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
/// \file EventAction.cc
/// \brief Implementation of the B2::EventAction class
#include "G4AnalysisManager.hh"
#include "EventAction.hh"
#include "StripHit.hh"
#include "StripSD.hh"

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{

  G4int id = event->GetEventID();
 // G4cout << "Event " << id+1 << "\r" << std::flush;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  auto hc = event->GetHCofThisEvent()->GetHC(0);
  G4int HitNb = hc -> GetSize();
  if (HitNb>0){
  G4int FireNb = 0;
  G4double AllEnergy[5][31] = {0};
  G4double FireTime = -1;
  G4double eTemp = 0;
  G4int Track = -1;
  for (int i = 0;i<HitNb;i++)
  { 
      StripHit* Hit = static_cast<StripHit*>(hc -> GetHit(i));
      FireNb = Hit -> GetStripNb();
      eTemp = Hit -> GetEdep();
      Track = Hit -> GetTrackID();
      if(Track>0&& Track<5) AllEnergy[Track-1][FireNb] = AllEnergy[Track-1][FireNb] + eTemp;
      else AllEnergy[4][FireNb] = AllEnergy[4][FireNb] + eTemp;
  }
  auto analysisManager = G4AnalysisManager::Instance();
  for(int i=0;i<5;i++)
  for(int j=0;j<31;j++)
  if (AllEnergy[i][j]>0){
  analysisManager -> FillNtupleDColumn(0,0,AllEnergy[i][j]/MeV);
  analysisManager -> FillNtupleDColumn(0,1,j);
  analysisManager -> FillNtupleDColumn(0,2,i+1);
  analysisManager -> AddNtupleRow(0);
  }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


