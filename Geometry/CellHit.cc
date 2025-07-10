//Custom Hit class taken from G4 Example basic/B4/B4c
#include "CellHit.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
// #include <cmath>
// #include <numbers>

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

void CellHit::SetTrackLength(G4double Trklen)
{
  fTrackLength = Trklen;
}

void CellHit::SetProcName(G4String procName)
{
  fProcName = procName;
}

void CellHit::SetProcId(G4int procId)
{
  fProcId = procId;
}

void CellHit::SetScatAngle(G4double theta)
{
  fScat = theta;
}

void CellHit::SetEta(G4double eta)
{
  fEta = eta;
}

void CellHit::SetPolarisation(G4ThreeVector pol)
{
  fPol = pol;
}


void CellHit::Print()
{
  G4cout << "Edep: " << std::setw(7) << G4BestUnit(fEdep, "Energy") << "|"
         << " track length: " << std::setw(7) << G4BestUnit(fTrackLength, "Length") << "|"
	 << " position: " << std::setw(7) << G4BestUnit(fPosition, "Length") << "|"
	 << " time: " << std::setw(7) << G4BestUnit(fTime, "Time") << "|"
	 << " DetId: " << std::setw(7) << detId << "|"
	 << " ProcessName: " << std::setw(7) << fProcName << "|"
	 << " ScatAngle: " << std::setw(7) << fScat*(180/CLHEP::pi) << "|"
	 << " eta: " << std::setw(7) << fEta*(180/CLHEP::pi) << G4endl;
    //<< "eta: " << std::setw(7) << fEta*(180/CLHEP::pi) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
