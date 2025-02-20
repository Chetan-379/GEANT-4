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
  
  G4cout << "=====================event No.: " << ievent << "====================" << G4endl;

}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  G4cout << "no. of hits: " << fX << G4endl;
  G4cout << "Energy deposition: " << fEdep << "\n" << G4endl;
  
  G4AnalysisManager *man = G4AnalysisManager::Instance(); 

  man->FillNtupleDColumn(1, 0, fX);
  man->AddNtupleRow(1);

  man->FillNtupleDColumn(3, 0, fEdep);
  man->AddNtupleRow(3);

  for (G4int i=0; i< Compt_edep.size(); i++){
    // G4cout << "Energy deposited in compton scattering: " << Compt_edep[i] << G4endl;
    compt_total_edep += Compt_edep[i];
  }
  
  G4cout << "total energy deposited via compton: " << compt_total_edep << G4endl;
  
  for (G4int i=0; i< Photo_edep.size(); i++){
    G4cout << "Energy deposited in PhotoElectric effect: " << Photo_edep[i] << G4endl;
  }

  G4double edep_ComptPhot = fEdep - compt_total_edep - Photo_edep[0];
  G4cout << "difference between fEdep and com+pho_edep: " << edep_ComptPhot << G4endl;


  runObject->X = Xarray;
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  runObject->Total_E = fEdep;
  runObject->Total_Compt_Edep = compt_total_edep;
  runObject->Photo_Edep = Photo_edep[0];
  runObject->diff_edep_ComptPhoto = edep_ComptPhot;
  runObject->tree->Fill();

  
  // G4cout << "size of Xarray is: " << Xarray.size() << G4endl;
  // G4cout << "size of Yarray is: " << Yarray.size() << G4endl;
  // G4cout << "size of Zarray is: " << Zarray.size() << G4endl;
  // G4cout << "size of Earray is: " << Earray.size() << G4endl;
  
  
  // G4cout << "elements in Xarray: " << G4endl;
  // for (G4int i=0; i<Xarray.size(); i++){
  // G4cout << Xarray[i] << G4endl;
  // }
  // G4cout << G4endl;

  // G4cout << "elements in Yarray: " << G4endl;
  // for (G4int i=0; i<Yarray.size(); i++){
  // G4cout << Yarray[i] << G4endl;
  // }
  // G4cout << G4endl;

  // G4cout << "elements in Zarray: " << G4endl;
  // for (G4int i=0; i<Zarray.size(); i++){
  // G4cout << Zarray[i] << G4endl;
  // }
  // G4cout << G4endl;

  // G4cout << "elements in Energy array: " << G4endl;
  // for (G4int i=0; i<Earray.size(); i++){
  // G4cout << Earray[i] << G4endl;
  // }
  // G4cout << G4endl;

  
  ievent++;
}
