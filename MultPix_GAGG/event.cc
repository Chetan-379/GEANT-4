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
  fX = 0.;
  fEdep = 0.;
  Xarray.clear();
  Yarray.clear();
  Zarray.clear();
  Earray.clear();
  Compt_edep.clear();
  Photo_edep.clear();
  compt_total_edep=0;
  photo_total_edep=0;
  nCompt_Before_Photo = 0;

  HitEdep_vec.clear();
  HitPosX_vec.clear();
  HitPosY_vec.clear();
  HitPosZ_vec.clear();
  HitTime_vec.clear();
  HitDetId_vec.clear();
  HitGunId_vec.clear();
  HitScatAngle_vec.clear();
  HitEta_vec.clear();
  HitEin_vec.clear();
  HitEout_vec.clear();
  HitProcId_vec.clear();

  for (int i = 0; i<8; i++){
    for (int j = 0; j<8; j++){
      matrix1[i][j] = 0;
      matrix2[i][j] = 0;
      
      matrix1_gen[i][j] = 0;
      matrix2_gen[i][j] = 0;

      runObject->Module1[i][j][nOpPho]=0;
      runObject->Module2[i][j][nOpPho]=0;

      runObject->Module1[i][j][nGenOp]=0;
      runObject->Module2[i][j][nGenOp]=0;

      runObject->Module1[i][j][Edep]=0;
      runObject->Module2[i][j][Edep]=0;

      runObject->Module1[i][j][E_elec]=0;
      runObject->Module2[i][j][E_elec]=0;
    }   
  }

   runObject->Edep_truth = 0;
   runObject->nPtcl_out = 0;
   runObject->E_outPtcl =0;
   runObject->E_verify =0;
   runObject->nSec_e = 0;
   runObject->nSec_pho = 0;
   
   runObject->E_eOut.clear();
   runObject->E_phoOut.clear();

   G4cout << "\n=====================event No.: " << ievent << "====================" << G4endl;    
   // for (int i =0; i<tgt_evt.size(); i++){
   //   if (ievent == tgt_evt[i]) {     
   //     G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  
   //   }     
   //   //else  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
   // }
   
   
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  if (Photo_edep.size() < 1) Photo_edep.push_back(0);
  
  for (G4int i=0; i< Compt_edep.size(); i++){
    compt_total_edep += Compt_edep[i];
  }

  for (G4int i=0; i< Photo_edep.size(); i++){
  }
  

  //****************Processing Hit collection of the Event*******************************
  // Get hits collections IDs (only once)
  
  if (fScintHCID == -1) {
    fScintHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintHitsCollection");
  }
  
  auto EvtHC = event->GetHCofThisEvent()->GetHC(0);

  
  G4int nHits = 0;

  for (int i =0; i< EvtHC->GetSize(); i++){
    CellHit* ScintHit =  dynamic_cast<CellHit*>(EvtHC->GetHit(i));
    nHits++;
      
    HitEdep_vec.push_back(ScintHit->GetEdep());
    HitPosX_vec.push_back((ScintHit->GetPosition())[0]);
    HitPosY_vec.push_back((ScintHit->GetPosition())[1]);
    HitPosZ_vec.push_back((ScintHit->GetPosition())[2]);
    HitTime_vec.push_back(ScintHit->GetTime());
    HitDetId_vec.push_back(ScintHit->GetDetID());
    HitGunId_vec.push_back(ScintHit->GetGunID());
    HitProcId_vec.push_back(ScintHit->GetProcId());

    HitScatAngle_vec.push_back(ScintHit->GetScatAngle());

    HitEta_vec.push_back(ScintHit->GetEta());

    HitEin_vec.push_back(ScintHit->GetEin());
    HitEout_vec.push_back(ScintHit->GetEout());

  }

  if (ievent == chkEvt) G4EventManager::GetEventManager()->KeepTheCurrentEvent();    
  
  //**************************************************************************************************

  //storing informations in the Tree
  runObject->X = Xarray;
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  //runObject->Total_E = fEdep;
  runObject->Total_E = compt_total_edep + Photo_edep[0];
  runObject->Total_Compt_Edep = compt_total_edep;
  runObject->Photo_Edep = Photo_edep[0];

  runObject->HitEdep = HitEdep_vec;
  runObject->HitPosX = HitPosX_vec;
  runObject->HitPosY = HitPosY_vec;
  runObject->HitPosZ = HitPosZ_vec;
  runObject->HitTime = HitTime_vec;
  runObject->HitDetId = HitDetId_vec;
  runObject->HitGunId = HitGunId_vec;
  runObject->HitProcId = HitProcId_vec;
  runObject->HitScatAngle = HitScatAngle_vec;
  runObject->HitEta = HitEta_vec;
  runObject->HitEout = HitEout_vec;
  runObject->HitEin = HitEin_vec;

  runObject->E_verify = runObject->E_outPtcl + runObject->Edep_truth;

  // runObject->E_elec_out = E_eOut_vec;
  // runObject->E_pho_out = E_phoOut_vec;
  
  runObject->tree->Fill();

  //G4cout << "\n\nsize of the matrix is: " << sizeof(matrix1) / sizeof(matrix1[0][0]) << G4endl;
  // G4cout << "\v============= right module Edep: =============" << G4endl;
  //  for (int i = 0; i < 8; i++) {          // Loop over rows
  //       for (int j = 0; j < 8; j++) {      // Loop over columns
  // 	  G4int precision = 0;
  // 	  if (runObject->Module1[i][j][Edep] ==0) precision = 0;
  // 	  if (runObject->Module1[i][j][Edep] !=0) precision = 3;
  // 	  G4cout << std::fixed << std::setprecision(precision) << std::setw(6) << runObject->Module1[i][j][Edep] << " ";
  //       }
  //       G4cout << G4endl;                      // New line after each row
  //   }

  //    G4cout << "\n\n\v============= right module nOpPho: =============" << G4endl;
  //  for (int i = 0; i < 8; i++) {          // Loop over rows
  //       for (int j = 0; j < 8; j++) {      // Loop over columns
  // 	  G4int precision = 0;
  // 	  if (runObject->Module1[i][j][nOpPho] ==0) precision = 0;
  // 	  if (runObject->Module1[i][j][nOpPho] !=0) precision = 1;
  // 	  G4cout << std::fixed << std::setprecision(precision) << std::setw(6) << runObject->Module1[i][j][nOpPho] << " ";
  //       }
  //       G4cout << G4endl;                      // New line after each row
  //   }

  // G4cout << "\n\n\v============= right module  elec_E: =============" << G4endl;
  //  for (int i = 0; i < 8; i++) {          // Loop over rows
  //       for (int j = 0; j < 8; j++) {      // Loop over columns
  // 	  G4int precision = 0;
  // 	  if (runObject->Module1[i][j][E_elec] ==0) precision = 0;
  // 	  if (runObject->Module1[i][j][E_elec] !=0) precision = 3;
  // 	  G4cout << std::fixed << std::setprecision(precision) << std::setw(6) << runObject->Module1[i][j][E_elec] << " ";
  //       }
  //       G4cout << G4endl;                      // New line after each row
  //   }
  
  //if (ievent == chkEvt) G4EventManager::GetEventManager()->KeepTheCurrentEvent();
  // for (int i =0; i<tgt_evt.size(); i++){
  //   if (ievent == tgt_evt[i]) {     
  //     G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");     
  //   }        
  // }
  ievent++;
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;  
}
