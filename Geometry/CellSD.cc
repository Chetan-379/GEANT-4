#include "CellSD.hh"

// #include "G4HCofThisEvent.hh"
// #include "G4SDManager.hh"
// #include "G4Step.hh"

// //namespace B4c
// //{

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

  // Create hits
  // fNofCells for cells + one more for total sums
  //for (G4int i = 0; i < fNofCells + 1; i++) {
  //fHitsCollection->insert(new CellHit());
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CellSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  //G4cout << "Energy Deposited in the step: " << edep << G4endl;
  // step length
  G4double stepLength = 0.;
  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
    stepLength = step->GetStepLength();
  }

  if (edep == 0. && stepLength == 0.) return false;

  auto touchable = (step->GetPreStepPoint()->GetTouchable());
  //G4int layerNumber = 0;
  //G4cout << "detector pos: " << touchable->GetVolume()->GetTranslation() << G4endl;

  // Get calorimeter cell id
  //auto layerNumber = touchable->GetReplicaNumber(1);
  //auto layerNumber = touchable->GetCopyNumber(1);
  // if(touchable->GetVolume() != nullptr) {layerNumber = touchable->GetVolume()->GetCopyNo();
  //   G4cout << "\n\nlayerNo.: " << layerNumber << G4endl;}
  //auto layerNumber =0;
  // auto test_var = step->GetTrack()->GetVolume()->GetCopyNo();

  // G4cout << "\ntestVar: " << test_var << G4endl;

   //getting the position of compton scattering point
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4ThreeVector compt_scat_point;
  G4int ParentID = step->GetTrack()->GetParentID();
  
  // G4int nCompt = 0;
  // if(proc == "compt" && ParentID == 0) nCompt++;
  // G4cout << "\n nCompt: " << nCompt << G4endl;
    
  
  // Get hit accounting data for this cell
  
  //auto hit =  (*fHitsCollection)[nCompt]; //filling the 0th element of hit collection
  
  auto hit = new CellHit();
  // if (!hit) {
  //   G4ExceptionDescription msg;
  //   msg << "Cannot access hit ";
  //   G4Exception("CellSD::ProcessHits()", "MyCode0004", FatalException, msg);
  // }

  // Get hit for total accounting
  //auto hitTotal = (*fHitsCollection)[fHitsCollection->entries() - 1];

  // Add values
  hit->Add(edep, stepLength);
  //hitTotal->Add(edep, stepLength);

 
  //set position and time
  G4ThreeVector DetPos = touchable->GetVolume()->GetTranslation();
  G4double currentTime = step->GetPreStepPoint()->GetGlobalTime();
  G4double TrkLen = step->GetTrack()->GetTrackLength();
  //G4double currentTime = step->GetPostStepPoint()->GetGlobalTime();
  //G4double currentTime = step->GetTrack()->GetGlobalTime();
  G4int DetId = touchable->GetVolume()->GetCopyNo();
  
  if (proc == "compt" && ParentID == 0) {
    compt_scat_point = step->GetPostStepPoint()->GetPosition();
    G4cout << "Compton scattering point: " << compt_scat_point << G4endl;
    //hit->SetPosition(DetPos);
    hit->SetPosition(compt_scat_point);
    hit->SetTime(currentTime);
    hit->SetDetectorID(DetId);
    hit->SetTrackLength(TrkLen);

    fHitsCollection->insert(hit);
  }

  //if(hit->GetTime()>0) fHitsCollection->insert(hit);
  //G4int id = fHitsCollection->insert(hit);
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
  //if(true){
    auto nofHits = fHitsCollection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << nofHits
           << " hits in the tracker chambers: " << G4endl;
    for (std::size_t i = 0; i < nofHits; ++i)
      (*fHitsCollection)[i]->Print();
  }
  
  nCompt = 0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//}  // namespace B4c
