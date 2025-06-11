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
  G4cout << "   Cell: total energy: " << std::setw(7) << G4BestUnit(CellEdep, "Energy")
         << "       total track length: " << std::setw(7) << G4BestUnit(CellTrackLength, "Length")
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
  nOptPho = 0;
  OptPho_PosX.clear();
  OptPho_PosY.clear();
  OptPho_PosZ.clear();
  OptPho_Energy.clear();
  OptPho_time.clear();
  
  TrkOnDet = 0;
  
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
  runObject->tree->Fill();  

  for (G4int i =0; i < OptPho_PosZ.size(); i++){
    if (OptPho_PosZ[i] < -1)  {G4EventManager::GetEventManager()->KeepTheCurrentEvent();
    }
  }


  if (fCellHCID == -1) {
    fCellHCID = G4SDManager::GetSDMpointer()->GetCollectionID("CellHitsCollection");
    
    //fGapHCID = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }

  //cout << 
  
  // Get hits collections
  //G4cout << "fCellHCID: " << fCellHCID << G4endl;
  auto CellHC = GetHitsCollection(G4int(fCellHCID), event);
  
  //auto gapHC = GetHitsCollection(fGapHCID, event);
  
  // Get hit with total values
  auto CellHit = (*CellHC)[CellHC->entries() - 1];
  //auto gapHit = (*gapHC)[gapHC->entries() - 1];

  
  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ((printModulo > 0) && (eventID % printModulo == 0)) {
    PrintEventStatistics(CellHit->GetEdep(), CellHit->GetTrackLength());
    G4cout << "--> End of event: " << eventID << "\n" << G4endl;
    //}

  
  
  
  ievent++;
  // G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;  
}
