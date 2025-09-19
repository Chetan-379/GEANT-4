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
  
  
  //G4cout << "=====================event No.: " << ievent << "====================" << G4endl;
 
  //if (ievent == chkEvt) G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  //else G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");
}

void MyEventAction::EndOfEventAction(const G4Event* event)
{
  if (Photo_edep.size() < 1) Photo_edep.push_back(0);
  
  for (G4int i=0; i< Compt_edep.size(); i++){
    compt_total_edep += Compt_edep[i];
  }  


  if (ievent == chkEvt) G4EventManager::GetEventManager()->KeepTheCurrentEvent();    
     //}
  
  //**************************************************************************************************
    
  
  //G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;

  double progress = 10.0 * ievent / (1.0 * 10000);
  int k = int (progress);
  if (k > decade)
    G4cout << 10 * k << " %" << G4endl;
  decade = k;

  ievent++;
  
}
