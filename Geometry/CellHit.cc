#include "CellHit.hh"

#include "G4UnitsTable.hh"

#include <iomanip>

//namespace B4c
//{

G4ThreadLocal G4Allocator<CellHit>* CellHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4bool CellHit::operator==(const CellHit& right) const
{
  return (this == &right) ? true : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellHit::SetPosition(G4ThreeVector xyz)
{
  fPosition = xyz;
}

void CellHit::SetTime(G4double time)
{
  fTime = time;
}

void CellHit::SetDetectorID(G4int DetID)
{
  detId = DetID;
}


void CellHit::Print()
{
  G4cout << "Edep: " << std::setw(7) << G4BestUnit(fEdep, "Energy")
         << " track length: " << std::setw(7) << G4BestUnit(fTrackLength, "Length")
	 << " position: " << std::setw(7) << G4BestUnit(fPosition, "Length")
	 << " time: " << std::setw(7) << G4BestUnit(fTime, "Time")
	 << " DetId: " << std::setw(7) << detId << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//}  // namespace B4c
