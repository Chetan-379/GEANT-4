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
   runObject->nOpGAGG_QE =0;
   runObject->nOpLYSO_QE =0;
   
   runObject->Edep_G_truth =0;
   runObject->Edep_L_truth =0;

   for (int iqud = 0; iqud<4; iqud ++){
     runObject->nOpG_quad[iqud] = 0;
     runObject->nOpL_quad[iqud] = 0;
   }

   storeHit = true;
   runObject->Det_clr_Idx.clear();
   runObject->Det_col_Idx.clear();
   runObject->Det_row_Idx.clear();
   runObject->Edep_truth_vec.clear();   
 
   //G4cout << "\n=====================event No.: " << ievent << "====================" << G4endl;     
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  // for (int i =0; i< store_log.size(); i++){
  //   //fRunAction->col_idx.push_back(100*clrIdx)
  //}
  
  runObject->tree->Fill();  
  ievent++;
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;

  // for (int i =0; i< runObject->Det_clr_Idx.size(); i++){
  //   G4cout << "color: " << runObject->Det_clr_Idx[i] << G4endl;
  //   G4cout << "row: " << runObject->Det_row_Idx[i] << G4endl;
  //   G4cout << "column: " << runObject->Det_col_Idx[i] << G4endl;
  //   G4cout << "edep: " << runObject->Edep_truth_vec[i] << "\n\n";

  //   if (runObject->Edep_truth_vec[i]>0.513) G4cout << "Event having larger edep: " << ievent << G4endl;
  //   break;
  // }

  //G4cout << "Hit vec size: " << runObject->Det_clr_Idx.size() << G4endl;
}
