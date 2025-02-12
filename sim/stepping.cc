#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
  fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
  //G4StepPoint *preStepPoint = step->GetPreStepPoint();
  G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
  const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  
  // G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
  // if(volume != fScoringVolume)
  // return; 
  

  G4ThreeVector position = detectorConstruction->sensDet->posDetector;
  G4double Energy_dep = detectorConstruction->sensDet->edep;
  G4double Xcord = position[0];
  G4double Ycord = position[1];
  G4double Zcord = position[2];

  //G4double edep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edep);

  fEventAction->AddX(Xcord);
  fEventAction->StorePos(Xcord, Ycord, Zcord, Energy_dep);

}
