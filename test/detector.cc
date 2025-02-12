#include "detector.hh"

MySensitiveDetector::MySensitiveDetector(G4String name) :
  G4VSensitiveDetector(name)
{}

MySensitiveDetector::~MySensitiveDetector()
{}


G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
  G4Track *track = aStep->GetTrack();
  track->SetTrackStatus(fStopAndKill);

  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  //G4ThreeVector posPhoton = preStepPoint->GetPosition();
  edep = aStep->GetTotalEnergyDeposit();

  //G4cout << "Photon position: " << posPhoton << G4endl;
  const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();

  G4int copyNo = touchable->GetCopyNumber();

  //G4cout << "Copy number: " << copyNo << G4endl;

  G4VPhysicalVolume *physVol = touchable->GetVolume();
  posDetector = physVol->GetTranslation();

  G4cout << "Detector position: " << posDetector << G4endl;
  // random++;
  // G4cout << "hit No." << random << G4endl;
  
  //test_pos.push_back(posDetector);

  G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(0, evt);
  man->FillNtupleDColumn(1, posDetector[0]);
  man->FillNtupleDColumn(2, posDetector[1]);
  man->FillNtupleDColumn(3, posDetector[2]);

  man->AddNtupleRow(0);
  
  return 0;
}
