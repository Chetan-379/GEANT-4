#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{
  fX = 0.;
}

MyEventAction::~MyEventAction()
{}

 void MyEventAction::BeginOfEventAction(const G4Event*)
{
  fX = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  G4cout << "X coordinate: " << fX << G4endl;

  G4AnalysisManager *man = G4AnalysisManager::Instance();

  man->FillNtupleDColumn(1, 0, fX);

  man->AddNtupleRow(1);
}
