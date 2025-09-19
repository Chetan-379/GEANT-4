//taken from G4 Example basic/B4/B4c
#ifndef CELLSD_HH
#define CELLSD_HH

#include "CellHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"
#include "G4EventManager.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

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
  G4bool summary = false; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
