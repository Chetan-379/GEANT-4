#ifndef CELLSD_HH
#define CELLSD_HH

#include "CellHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "globals.hh"


class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//namespace cellSD
//{

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class CellSD : public G4VSensitiveDetector
{
  public:
    CellSD(const G4String& name, const G4String& hitsCollectionName, G4int nofCells);
    ~CellSD() override = default;

    // methods from base class
    void Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void EndOfEvent(G4HCofThisEvent* hitCollection) override;

  private:
    CellHitsCollection* fHitsCollection = nullptr;
    G4int fNofCells = 0;
  G4int nCompt = 0; 
};

//}  // namespace B4c

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
