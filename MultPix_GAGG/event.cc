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
  // G4cout << "   InStrip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
  //        << "       InStrip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length") << G4endl;
  // G4cout << "   Strip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
  //        << "   Strip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length")
  //        << G4endl;
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
  HitTrkLen_vec.clear();
  HitProcId_vec.clear();
  HitScatAngle_vec.clear();
  HitEta_vec.clear();
  HitPol0_vec.clear();
  HitPol1_vec.clear();
  HitPol2_vec.clear();
  HitEin_vec.clear();
  HitEout_vec.clear();
  HitScatMomX_vec.clear();
  HitScatMomY_vec.clear();
  HitScatMomZ_vec.clear();
  
  
  G4cout << "=====================event No.: " << ievent << "====================" << G4endl;
 
  //if (ievent == chkEvt) G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  //else G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  if (Photo_edep.size() < 1) Photo_edep.push_back(0);
  
  for (G4int i=0; i< Compt_edep.size(); i++){
    compt_total_edep += Compt_edep[i];
  }

  for (G4int i=0; i< Photo_edep.size(); i++){
  }
  

  // //****************Processing Hit collection of the Event*******************************

  // // // Get hits collections IDs (only once)
  // if (fInStripHCID == -1) {
  //   fInStripHCID = G4SDManager::GetSDMpointer()->GetCollectionID("InStripHitsCollection");
  //   fOutStripHCID = G4SDManager::GetSDMpointer()->GetCollectionID("StripHitsCollection");
  //   fOuterMostStripHCID = G4SDManager::GetSDMpointer()->GetCollectionID("OuterMostStripHitsCollection");
  // }
  
  // // Get hits collections
  // auto InStripHC = event->GetHCofThisEvent()->GetHC(0);
  // auto OutStripHC = event->GetHCofThisEvent()->GetHC(1);
  // auto OuterMostStripHC = event->GetHCofThisEvent()->GetHC(2);

  // auto ScintHC = new G4THitsCollection<CellHit>;

  // ///////////////////////////////
  // //Merging the Hit collections from innermost layer and middle 2layers and the outermost layer in a single hit collection  

  // for (int iHit = 0; iHit < InStripHC->GetSize(); iHit++){
  //   auto hit = dynamic_cast<CellHit*>(InStripHC->GetHit(iHit));
  //   ScintHC->insert(hit);
  //   //auto scatAngl = hit->GetScatAngle()*180/CLHEP::pi;
  //   //if(scatAngl > 89 && scatAngl < 91 && hit->GetDetID() == 1051) G4cout << "\n\n=====event with 90: " << ievent << G4endl;
    
  // }

  // for (int iHit = 0; iHit < OutStripHC->GetSize(); iHit++){
  //   auto hit = dynamic_cast<CellHit*>(OutStripHC->GetHit(iHit));
  //   ScintHC->insert(hit); 
  // }

  // for (int iHit = 0; iHit < OuterMostStripHC->GetSize(); iHit++){
  //   auto hit = dynamic_cast<CellHit*>(OuterMostStripHC->GetHit(iHit));
  //   ScintHC->insert(hit); 
  // }
  // ///////////////////////////////

  // //G4cout << "\ntotal size of the hit collection in this event: " << ScintHC->GetSize() << G4endl;
  
  // G4int nHits = 0;

  // for (int i =0; i< ScintHC->GetSize(); i++){
  //   CellHit* ScintHit =  dynamic_cast<CellHit*>(ScintHC->GetHit(i));
  //   nHits++;
      
  //   HitEdep_vec.push_back(ScintHit->GetEdep());
  //   HitPosX_vec.push_back((ScintHit->GetPosition())[0]);
  //   HitPosY_vec.push_back((ScintHit->GetPosition())[1]);
  //   HitPosZ_vec.push_back((ScintHit->GetPosition())[2]);
  //   HitTime_vec.push_back(ScintHit->GetTime());
  //   HitDetId_vec.push_back(ScintHit->GetDetID());
  //   HitGunId_vec.push_back(ScintHit->GetGunID());
  //   HitTrkLen_vec.push_back(ScintHit->GetTrackLength());
  //   HitProcId_vec.push_back(ScintHit->GetProcId());

  //   HitScatAngle_vec.push_back(ScintHit->GetScatAngle());

  //   HitEta_vec.push_back(ScintHit->GetEta());
  //   HitPol0_vec.push_back((ScintHit->GetPolarisation())[0]);
  //   HitPol1_vec.push_back((ScintHit->GetPolarisation())[1]);
  //   HitPol2_vec.push_back((ScintHit->GetPolarisation())[2]);

  //   HitEin_vec.push_back(ScintHit->GetEin());
  //   HitEout_vec.push_back(ScintHit->GetEout());

  //   HitScatMomX_vec.push_back((ScintHit->GetScatMom())[0]);
  //   HitScatMomY_vec.push_back((ScintHit->GetScatMom())[1]);
  //   HitScatMomZ_vec.push_back((ScintHit->GetScatMom())[2]);

  // }

  // if (ievent == chkEvt) G4EventManager::GetEventManager()->KeepTheCurrentEvent();    
  //    //}
  
  //**************************************************************************************************

  //storing informations in the Tree
  runObject->X = Xarray;
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  runObject->Total_E = fEdep;
  runObject->Total_Compt_Edep = compt_total_edep;
  runObject->Photo_Edep = Photo_edep[0];

  // runObject->HitEdep = HitEdep_vec;
  // runObject->HitPosX = HitPosX_vec;
  // runObject->HitPosY = HitPosY_vec;
  // runObject->HitPosZ = HitPosZ_vec;
  // runObject->HitTime = HitTime_vec;
  // runObject->HitTrklen = HitTrkLen_vec;
  // runObject->HitDetId = HitDetId_vec;
  // runObject->HitGunId = HitGunId_vec;
  // runObject->HitProcId = HitProcId_vec;
  // runObject->HitScatAngle = HitScatAngle_vec;
  // runObject->HitEta = HitEta_vec;
  // runObject->HitPol0 = HitPol0_vec;
  // runObject->HitPol1 = HitPol1_vec;
  // runObject->HitPol2 = HitPol2_vec;
  // runObject->HitEout = HitEout_vec;
  // runObject->HitEin = HitEin_vec;
  // runObject->HitScatMomX = HitScatMomX_vec;
  // runObject->HitScatMomY = HitScatMomY_vec;
  // runObject->HitScatMomZ = HitScatMomZ_vec;
  
  
  runObject->tree->Fill();
    
  ievent++;
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;  
}
