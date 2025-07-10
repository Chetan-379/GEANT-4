//costum hit class taken from G4 Example basic/B4/B4c
#ifndef CellHit_h
#define CellHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "globals.hh"
#include "CLHEP/Units/PhysicalConstants.h"

class CellHit : public G4VHit
{
public:
  CellHit() = default;
  CellHit(const CellHit&) = default;
  ~CellHit() override = default;
  
  // operators
  CellHit& operator=(const CellHit&) = default;
  G4bool operator==(const CellHit&) const;
  
  inline void* operator new(size_t);
  inline void operator delete(void*);
  
  // methods from base class
  void Draw() override {}
  void Print() override;
  
  // methods to handle data
  void Add(G4double de, G4double dl);
  void SetPosition(G4ThreeVector xyz);
  void SetTime(G4double time);
  void SetDetectorID(G4int DetID);
  void SetTrackLength(G4double Trklen);
  void SetProcName(G4String procName);
  void SetProcId(G4int procId);
  void SetScatAngle(G4double theta);
  void SetEta(G4double eta);
  void SetPolarisation(G4ThreeVector pol);
    
  
  // get methods
  G4double GetEdep() const;
  G4double GetTrackLength() const;
  G4ThreeVector GetPosition() const;
  G4double GetTime() const;
  G4int GetDetID() const;
  G4String GetProcName() const;
  G4int GetProcId() const;
  G4double GetScatAngle() const;
  G4double GetEta() const;
  G4ThreeVector GetPolarisation() const;

  
private:
  G4double fEdep = 0.;  ///< Energy deposit in the sensitive volume
  G4double fTrackLength = 0.;  ///< Track length in the  sensitive volume
  G4ThreeVector fPosition, fPol;
  G4double fTime = 0;
  G4int detId = -1;
  G4String fProcName = " ";
  G4int fProcId = -1;
  G4double fScat = 10000.;
  G4double fEta = 10000.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using CellHitsCollection = G4THitsCollection<CellHit>;

extern G4ThreadLocal G4Allocator<CellHit>* CellHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CellHit::operator new(size_t)
{
  if (!CellHitAllocator) {
    CellHitAllocator = new G4Allocator<CellHit>;
  }
  void* hit;
  hit = (void*)CellHitAllocator->MallocSingle();
  return hit;
}

inline void CellHit::operator delete(void* hit)
{
  if (!CellHitAllocator) {
    CellHitAllocator = new G4Allocator<CellHit>;
  }
  CellHitAllocator->FreeSingle((CellHit*)hit);
}

inline void CellHit::Add(G4double de, G4double dl)
{
  fEdep += de;
  fTrackLength += dl;
}

inline G4double CellHit::GetEdep() const
{
  return fEdep;
}

inline G4double CellHit::GetTrackLength() const
{
  return fTrackLength;
}

inline G4ThreeVector CellHit::GetPosition() const
{
  return fPosition;
}

inline G4double CellHit::GetTime() const
{
  return fTime;
}

inline G4int CellHit::GetDetID() const
{
  return detId;
}

inline G4String CellHit::GetProcName() const
{
  return fProcName;
}

inline G4int CellHit::GetProcId() const
{
  return fProcId;
}

inline G4double CellHit::GetScatAngle() const
{
  return fScat;
}

inline G4double CellHit::GetEta() const
{
  return fEta;
}

inline G4ThreeVector CellHit::GetPolarisation() const
{
  return fPol;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
