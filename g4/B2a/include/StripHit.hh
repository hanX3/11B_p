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
/// \file TrackerHit.hh
/// \brief Definition of the B2::TrackerHit class

#ifndef B2StripHit_h
#define B2StripHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"



/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class StripHit : public G4VHit
{
  public:
    StripHit();
    StripHit(const StripHit&) = default;
    ~StripHit() override;

    // operators
    StripHit& operator=(const StripHit&) = default;
    G4bool operator==(const StripHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
 //   void Draw() override;
  //  void Print() override;

    // Set methods
    void SetTrackID  (G4int track)      { fTrackID = track; };
    void SetStripNb(G4int strip)      { fStripNb = strip; };
    void SetEdep     (G4double de)      { fEdep = de; };
    void SetTime   (G4double time)   {fTime = time;};

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetStripNb() const   { return fStripNb; };
    G4double GetEdep() const     { return fEdep; };
    G4double GetTime() const {return fTime;};
   // G4ThreeVector GetPos() const { return fPos; };

  private:
    G4int         fTrackID = -1;
    G4int         fStripNb = -1;
    G4double      fEdep = 0.;
    G4double      fTime = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<StripHit> StripHitsCollection;

extern G4ThreadLocal G4Allocator<StripHit>* StripHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* StripHit::operator new(size_t)
{
  if(!StripHitAllocator)
      StripHitAllocator = new G4Allocator<StripHit>;
  return (void *) StripHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void StripHit::operator delete(void *hit)
{
  StripHitAllocator->FreeSingle((StripHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



#endif
