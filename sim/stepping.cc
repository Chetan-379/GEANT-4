#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
  fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
  G4StepPoint *preStepPoint = step->GetPreStepPoint();
  G4ThreeVector position = preStepPoint->GetPosition();
  G4double Xcord = position[0];
  G4double Ycord = position[1];
  G4double Zcord = position[2];

  fEventAction->AddX(Xcord);
  //fEventAction->StoreX(Xcord);
  //fEventAction->FillTree();
  // fEventAction->StoreY(Ycord);
  // fEventAction->StoreZ(Zcord);
  //fEventAction->StoreHitPos(position);
}
