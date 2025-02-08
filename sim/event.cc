#include "event.hh"

MyEventAction::MyEventAction(MyRunAction* runAction)
{
  fX = 0.;
  //fill(HitsArray.begin(), HitsArray.end(), 0);
  //TTree *tree;
  //runObject.tree->Branch("positionX", &Xarray, "positionX/I");

  //MyRunAction runObject;
  //runObject.x= 5000;
  //int p=90;
  //runObject.FillTree();
  runObject = runAction;
  
}

MyEventAction::~MyEventAction()
{}

 void MyEventAction::BeginOfEventAction(const G4Event*)
{
  fX = 0.;
  //Xarray->clear();
  Xarray.clear();

  //MyRunAction *runObject;
  //runObject->tree;
}

void MyEventAction::EndOfEventAction(const G4Event*)
{
  G4cout << "no. of hits: " << fX << G4endl;
  //G4cout << "size of X array:  " << Xarray->size() << G4endl;
  // G4cout << "size of Y array:  " << Yarray.size() << G4endl;
  // G4cout << "size of Z array:  " << Zarray.size() << G4endl;
   //G4cout << "size of Hits array:  " << HitsArray.size() << G4endl;
  
  G4AnalysisManager *man = G4AnalysisManager::Instance(); 

  man->FillNtupleDColumn(1, 0, fX);
  man->AddNtupleRow(1);

  //runObject.tree->Fill();
  runObject->p = Xarray;
  runObject->tree->Fill();
  
  //runObject.FillTree(Xarray);
  //tree->Print();

  // TTree *tree;

  //runObject.tree = ;
  //runObject.Xcoord = Xarray;
  //runObject.tree->Branch("positionX", &Xarray, "positionX/I");
  //tree1->Fill();
  //  if (Xarray->size()>0){ 
  G4cout << "size of Xarray is: " << Xarray.size() << G4endl;
  G4cout << "X positions in event: " << ievent << G4endl;
  for (G4int i=0; i<Xarray.size(); i++){
  G4cout << Xarray[i] << G4endl;
    // runObject.Xcoord = Xarray[i];
    // runObject.tree->Fill()
  //MyRunAction runObject;
  //runObject.FillTree(Xarray);
  //}
  //}
    //man->FillNtupleDColumn(2, 0, Xarray[0]);
   ///   tree->Branch("Xcoord",&Xarray[i]);
  
  }
  //man->AddNtupleRow(2);
  G4cout << G4endl;
  ievent++;
}
