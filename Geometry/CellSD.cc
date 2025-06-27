//taken from G4 example basic/B4/B4c
#include "CellSD.hh"

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellSD::CellSD(const G4String& name, const G4String& hitsCollectionName,
	       G4int nofCells)
  : G4VSensitiveDetector(name), fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new CellHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CellSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;

  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
    stepLength = step->GetStepLength();
  }

  if (edep == 0. && stepLength == 0.) return false;

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  //getting the position of compton scattering point
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4ThreeVector compt_scat_point;
  G4int ParentID = step->GetTrack()->GetParentID();
  
  // G4int nCompt = 0;
  // if(proc == "compt" && ParentID == 0) nCompt++;
  // G4cout << "\n nCompt: " << nCompt << G4endl;
    
  
  // Get hit accounting data for this cell
  auto hit = new CellHit();

  // Add values
  hit->Add(edep, stepLength);
 
  //set position, time etc
  G4double currentTime = step->GetPreStepPoint()->GetGlobalTime();
  G4double TrkLen = step->GetTrack()->GetTrackLength();
  G4int DetId = touchable->GetVolume()->GetCopyNo();
  G4double scat_theta = -1;
  G4ThreeVector vin = step->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vout = step->GetPostStepPoint()->GetMomentumDirection();

  // G4ThreeVector vin = step->GetPreStepPoint()->GetMomentum();
  // G4ThreeVector vout = step->GetPostStepPoint()->GetMomentum();

  // G4cout << "vin mag" << vin.mag() << G4endl;
  if ((proc == "compt" || proc == "phot") && ParentID == 0) {
    compt_scat_point = step->GetPostStepPoint()->GetPosition();

    G4cout << "\nCompton Scat point: " << compt_scat_point << G4endl;
    G4cout << "pre step momDir: " << vin << G4endl;
    G4cout << "post step momDir: " << vout << "\n\n";

    if (proc == "compt") scat_theta = acos(vin.dot(vout));
    
    
    hit->SetPosition(compt_scat_point);
    hit->SetTime(currentTime);
    hit->SetDetectorID(DetId);
    hit->SetTrackLength(TrkLen);
    hit->SetProcName(proc);
    hit->SetScatAngle(scat_theta);

    fHitsCollection->insert(hit);
  }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    auto nofHits = fHitsCollection->entries();
    
    G4cout << "\n" << "-------->Hit summary of the SD "  << G4endl;
    for (std::size_t i = 0; i < nofHits; ++i)
      (*fHitsCollection)[i]->Print();
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
