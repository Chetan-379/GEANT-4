#include "event.hh"

MyEventAction::MyEventAction(MyRunAction* runAction)
{
  fX = 0.;
  fEdep = 0.;
  runObject = runAction;  
}

MyEventAction::~MyEventAction()
{}

 void MyEventAction::BeginOfEventAction(const G4Event*)
{
  G4cout << "\n\n\t\t============Event Started=============" << G4endl;
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
  n_Refraction = 0;
  n_Reflection = 0;
  n_TIR = 0;  
  TrkOnDet = 0;
  
  //G4cout << "=====================event No.: " << ievent << "====================" << G4endl;

}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  if (Photo_edep.size() < 1) Photo_edep.push_back(0);
  
  
  G4AnalysisManager *man = G4AnalysisManager::Instance(); 
  
  man->FillNtupleDColumn(1, 0, fX);
  man->AddNtupleRow(1);

  man->FillNtupleDColumn(3, 0, fEdep);
  man->AddNtupleRow(3);

  for (G4int i=0; i< Compt_edep.size(); i++){
    compt_total_edep += Compt_edep[i];
  }
  
  
  for (G4int i=0; i< Photo_edep.size(); i++){
  }

  
  runObject->X = Xarray;  
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  runObject->Total_E = fEdep;
  runObject->Total_Compt_Edep = compt_total_edep;
  runObject->Photo_Edep = Photo_edep[0];
  runObject->diff_edep_ComptPhoto = edep_ComptPhot;
  runObject->nOpticalPhotons = nOptPho;
  runObject->nOptPhoOnDetEnd = TrkOnDet;
  runObject->Opt_Photon_PosX = OptPho_PosX;
  runObject->Opt_Photon_PosY = OptPho_PosY;
  runObject->Opt_Photon_PosZ = OptPho_PosZ;
  runObject->Opt_Photon_Energy = OptPho_Energy;
  runObject->Opt_Photon_Time = OptPho_time;
  runObject->nRefraction = n_Refraction;
  runObject->nReflection = n_Reflection;
  runObject->nTIR = n_TIR;
  // auto EvtMan = G4EventManager::GetEventManager();
  // if (ievent == 9) {
  //   EvtMan->KeepTheCurrentEvent();
  // }

 //  if (ievent == 10) {
//   G4EventManager::GetEventManager()->KeepTheCurrentEvent();
// }

  for (G4int i =0; i < OptPho_PosZ.size(); i++){
    if (OptPho_PosZ[i] < -1)  {G4EventManager::GetEventManager()->KeepTheCurrentEvent();
      //G4cout << "\n\nOptPhoton_z: " << OptPho_PosZ[i] << G4endl;
      //G4cout << "Event no. is: " << ievent << G4endl;}
  }

  }

      


ievent++;
 
//G4cout << "\n\t\t\t****************END OF EVENT*******************\n" << G4endl;
 
 
}
