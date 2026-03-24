#include "event.hh"

MyEventAction::MyEventAction(MyRunAction* runAction)
{
  runObject = runAction;  
}

MyEventAction::~MyEventAction()
{}

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
  runObject->Det_nCompt.clear();
  runObject->Det_nPhoto.clear();
  runObject->Det_theta.clear();
   
  G4cout << "\n=====================event No.: " << ievent << "====================" << G4endl;
  if (ievent == chkEvt) G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");

  auto vtx1 = event->GetPrimaryVertex(0);
  auto mom1 = vtx1->GetPrimary(0)->GetMomentumDirection();

  runObject->h_momXvsY->Fill(mom1[0],mom1[1]);
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{ 
  runObject->tree->Fill();

  if (ievent == chkEvt) {
    G4EventManager::GetEventManager()->KeepTheCurrentEvent();
    G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
  }
  
  ievent++;

  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;

  for (int i =0; i< runObject->Det_clr_Idx.size(); i++){
    // G4cout << "color: " << runObject->Det_clr_Idx[i] << G4endl;
    // G4cout << "row: " << runObject->Det_row_Idx[i] << G4endl;
    // G4cout << "column: " << runObject->Det_col_Idx[i] << G4endl;
    // G4cout << "edep: " << runObject->Edep_truth_vec[i] << G4endl;
    // G4cout << "nCompt: " << runObject->Det_nCompt[i] << G4endl;
    // G4cout << "nPhoto: " << runObject->Det_nPhoto[i] << G4endl;
    // G4cout << "scat_theta: " << runObject->Det_theta[i] << "\n\n";
  }

  // G4cout << "scat theta vec size: " << runObject->Det_theta.size() << G4endl;
  // G4cout << "row Idx vec size: " << runObject->Det_row_Idx.size() << G4endl;
  
}
