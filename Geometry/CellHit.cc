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

void CellHit::Print()
{
  G4cout << "Edep: " << std::setw(7) << G4BestUnit(fEdep, "Energy")
         << " track length: " << std::setw(7) << G4BestUnit(fTrackLength, "Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//}  // namespace B4c
