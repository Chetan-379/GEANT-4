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
  // print event statistics
  G4cout << "   InStrip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
         << "       InStrip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length") << G4endl;
  G4cout << "   Strip: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
         << "   Strip total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length")
         << G4endl;

  //G4cout << "Printing Hits!!" << G4endl; 
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
  nOptPho = 0;
  OptPho_PosX.clear();
  OptPho_PosY.clear();
  OptPho_PosZ.clear();
  OptPho_Energy.clear();
  OptPho_time.clear();
  
  TrkOnDet = 0;

  InHitEdep_vec.clear();
  InHitPosX_vec.clear();
  InHitPosY_vec.clear();
  InHitPosZ_vec.clear();
  InHitTime_vec.clear();
  InHitDetId_vec.clear();
  InHitTrkLen_vec.clear();
  
  //G4cout << "=====================event No.: " << ievent << "====================" << G4endl;

}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  if (Photo_edep.size() < 1) Photo_edep.push_back(0);
  
  for (G4int i=0; i< Compt_edep.size(); i++){
    // G4cout << "Energy deposited in compton scattering: " << Compt_edep[i] << G4endl;
    compt_total_edep += Compt_edep[i];
  }

  for (G4int i=0; i< Photo_edep.size(); i++){
    //G4cout << "Energy deposited in PhotoElectric effect: " << Photo_edep[i] << G4endl;
  }
  
  runObject->X = Xarray;
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  runObject->Total_E = fEdep;
  runObject->Total_Compt_Edep = compt_total_edep;
  runObject->Photo_Edep = Photo_edep[0];
  //runObject->diff_edep_ComptPhoto = edep_ComptPhot;
  runObject->nOpticalPhotons = nOptPho;
  runObject->nOptPhoOnDetEnd = TrkOnDet;
  runObject->Opt_Photon_PosX = OptPho_PosX;
  runObject->Opt_Photon_PosY = OptPho_PosY;
  runObject->Opt_Photon_PosZ = OptPho_PosZ;
  runObject->Opt_Photon_Energy = OptPho_Energy;
  runObject->Opt_Photon_Time = OptPho_time;
  //runObject->tree->Fill();  

  for (G4int i =0; i < OptPho_PosZ.size(); i++){
    if (OptPho_PosZ[i] < -1)  {G4EventManager::GetEventManager()->KeepTheCurrentEvent();
    }
  }

  //========================================================making hit collection=================================================================================================================================
  // Get hits collections IDs (only once)
  if (fInStripHCID == -1) {
    fInStripHCID = G4SDManager::GetSDMpointer()->GetCollectionID("InStripHitsCollection");
    fStripHCID = G4SDManager::GetSDMpointer()->GetCollectionID("StripHitsCollection");   
  }

  //cout << 
  
  // Get hits collections
  //G4cout << "fCellHCID: " << fCellHCID << G4endl;
  auto InStripHC = GetHitsCollection(G4int(fInStripHCID), event);
  auto StripHC = GetHitsCollection(G4int(fStripHCID), event);

  //G4cout << "type of InStripHC: " << typeid(InStripHC).name() << G4endl;
  //G4cout << "type of StripHC: " << typeid(StripHC).name() << G4endl;
  //G4cout << "size CellHC: " << CellHC->entries() << G4endl;
  //G4cout << "size CellHC: " << event->GetHCofThisEvent()->GetNumberOfCollections() << G4endl;
  //G4cout << "size CellHC: " << event->GetHCofThisEvent()->GetHC(0)->GetSize() << G4endl;
  // std::vector<G4AttValue> *HitVars = event->GetHCofThisEvent()->GetHC(0)->GetHit(0)->CreateAttValues();
  // G4cout << "Total variables stored in the hit: " << HitVars->size() << G4endl;
  //G4cout << "Cell Hit: " << (event->GetHCofThisEvent()->GetHC(0)->GetHit(0)->CreateAttValues())[0].GetName() << G4endl;
  
  //auto gapHC = GetHitsCollection(fGapHCID, event);
  G4cout << "size of strip hits: " << StripHC->entries() << G4endl;
  for (int i =0; i< StripHC->entries(); i++){
    //auto test_Hit = InStripHC->GetHit(i);
    auto ScintHit = (*StripHC)[i];
    //auto Strip_Hit = (*StripHC) 
    //G4cout << "type of test_Hit is: " << typeid(test_Hit).name() << G4endl;
    //G4cout << "energy stored in Hit is: " << test_Hit->GetEdep() << G4endl;
    // G4ThreeVector position = ScintHit->GetPosition();
    // std::vector<double> position_stl = {position.x(), position.y(), position.z()};
    //ScintHit->Print();
    if(ScintHit->GetTrackLength() > 0 && ScintHit->GetTime() > 0){
      G4cout << "Edep: " << ScintHit->GetEdep() << G4endl;
      G4cout << "Time: " << ScintHit->GetTime() << G4endl;
      G4cout << "DetId: " << ScintHit->GetDetID() << G4endl;
      
      InHitEdep_vec.push_back(ScintHit->GetEdep());
      InHitPosX_vec.push_back((ScintHit->GetPosition())[0]);
      InHitPosY_vec.push_back((ScintHit->GetPosition())[1]);
      InHitPosZ_vec.push_back((ScintHit->GetPosition())[2]);
      InHitTime_vec.push_back(ScintHit->GetTime());
      InHitDetId_vec.push_back(ScintHit->GetDetID());
      InHitTrkLen_vec.push_back(ScintHit->GetTrackLength());
    }
  }
  
  // Get hit with total values
  auto InStripHit = (*InStripHC)[InStripHC->entries() - 1];
  auto StripHit = (*StripHC)[StripHC->entries() - 1];

  
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ((printModulo > 0) && (eventID % printModulo == 0)) {
    PrintEventStatistics(InStripHit->GetEdep(), InStripHit->GetTrackLength());
    PrintEventStatistics(StripHit->GetEdep(), StripHit->GetTrackLength());
    G4cout << "--> End of event: " << eventID << "\n" << G4endl;
    //}
    //==========================================================================================================================================================================================================


    runObject->InHitEdep = InHitEdep_vec;
    runObject->InHitPosX = InHitPosX_vec;
    runObject->InHitPosY = InHitPosY_vec;
    runObject->InHitPosZ = InHitPosZ_vec;
    runObject->InHitTime = InHitTime_vec;
    runObject->InHitTrklen = InHitTrkLen_vec;
    runObject->InHitDetId = InHitDetId_vec;
    runObject->tree->Fill();  
  
  
  ievent++;
  // G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;  
}
