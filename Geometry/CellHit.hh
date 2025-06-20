#ifndef B4cCalorHit_h
#define B4cCalorHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "globals.hh"

//namespace B4c
//{

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

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
    
  
  // get methods
  G4double GetEdep() const;
  G4double GetTrackLength() const;
  G4ThreeVector GetPosition() const;
  G4double GetTime() const;
  G4int GetDetID() const;
  

  
  private:
  G4double fEdep = 0.;  ///< Energy deposit in the sensitive volume
  G4double fTrackLength = 0.;  ///< Track length in the  sensitive volume
  G4ThreeVector fPosition;
  G4double fTime;
  G4int detId;
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

//inline void CellHit::SetPosition(G4ThreeVector xyz);

//inline void CellHit::SetDetectorID(G4int DetID);

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


//}  // namespace B4c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
