#include "event.hh"

MyEventAction::MyEventAction(MyRunAction* runAction)
{
  runObject = runAction;  
}

MyEventAction::~MyEventAction()
{}


void MyEventAction::BeginOfEventAction(const G4Event* event)
{
  
  //G4cout << "=====================event No.: " << ievent << "====================" << G4endl;
 
  //if (ievent == chkEvt) G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  //else G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");

  nOpPhotons = 0;
  nOpPhotons_end = 0;
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  if (ievent == chkEvt) G4EventManager::GetEventManager()->KeepTheCurrentEvent();

  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(2, 0, nOpPhotons);
  man->FillNtupleIColumn(2, 1, nOpPhotons_end);
  man->AddNtupleRow(2);
	  
    
  
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;

  
}
