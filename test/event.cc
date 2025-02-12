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

}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  G4cout << "no. of hits: " << fX << G4endl;
  G4cout << "Energy deposition: " << fEdep << G4endl;
  
  G4AnalysisManager *man = G4AnalysisManager::Instance(); 

  man->FillNtupleDColumn(1, 0, fX);
  man->AddNtupleRow(1);

  runObject->X = Xarray;
  runObject->Y = Yarray;
  runObject->Z = Zarray;
  runObject->E = Earray;
  runObject->Total_E = fEdep;
  runObject->tree->Fill();
  
  // G4cout << "size of Xarray is: " << Xarray.size() << G4endl;
  // G4cout << "size of Yarray is: " << Yarray.size() << G4endl;
  // G4cout << "size of Zarray is: " << Zarray.size() << G4endl;
  // G4cout << "size of Earray is: " << Earray.size() << G4endl;
  
  
  G4cout << "elements in Xarray: " << G4endl;
  for (G4int i=0; i<Xarray.size(); i++){
  G4cout << Xarray[i] << G4endl;
  }
  G4cout << G4endl;

  G4cout << "elements in Yarray: " << G4endl;
  for (G4int i=0; i<Yarray.size(); i++){
  G4cout << Yarray[i] << G4endl;
  }
  G4cout << G4endl;

  G4cout << "elements in Zarray: " << G4endl;
  for (G4int i=0; i<Zarray.size(); i++){
  G4cout << Zarray[i] << G4endl;
  }
  G4cout << G4endl;

  G4cout << "elements in Energy array: " << G4endl;
  for (G4int i=0; i<Earray.size(); i++){
  G4cout << Earray[i] << G4endl;
  }
  G4cout << G4endl;

  
  ievent++;
}
