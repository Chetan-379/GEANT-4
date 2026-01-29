#include "event.hh"

MyEventAction::MyEventAction(MyRunAction* runAction)
{
  fX = 0.;
  fEdep = 0.;
  runObject = runAction;  
}

MyEventAction::~MyEventAction()
{}

CellHitsCollection* MyEventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  auto hitsCollection = static_cast<CellHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  if (!hitsCollection) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

void MyEventAction::PrintEventStatistics(G4double CellEdep, G4double CellTrackLength) const
{
  G4cout << "   InStrip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
         << "       InStrip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length") << G4endl;
  G4cout << "   Strip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
         << "   Strip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length")
         << G4endl;
}

void MyEventAction::BeginOfEventAction(const G4Event* event)
{  
   runObject->Edep_truth = 0;
   runObject->nOpGAGG =0;
   runObject->nOpLYSO =0;
   runObject->nOpGAGG_truth =0;
   runObject->nOpLYSO_truth =0;
   runObject->Edep_G_truth =0;
   runObject->Edep_L_truth =0;
 
   G4cout << "\n=====================event No.: " << ievent << "====================" << G4endl;     
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  // G4cout << "\nGAGG_Edep: " << Edep_G_truth << G4endl;
  // G4cout << "\nLYSO_Edep: " << Edep_L_truth << G4endl;
  
  runObject->tree->Fill();
  ievent++;
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;  
}
