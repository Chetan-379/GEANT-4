#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
  fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
  G4StepPoint *preStepPoint = step->GetPreStepPoint();
  G4StepPoint *postStepPoint = step->GetPostStepPoint();
    
  G4ThreeVector phopos_pre = preStepPoint->GetPosition();
  G4ThreeVector phopos_post = postStepPoint->GetPosition();

  G4double edep = step->GetTotalEnergyDeposit();

  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  // G4String pre_proc;
  // if (step->GetPreStepPoint()->GetProcessDefinedStep()){
  //   pre_proc= step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName();
  // }
    
  G4double KE_i = preStepPoint->GetKineticEnergy();
  G4double KE_f = postStepPoint->GetKineticEnergy();

  G4ThreeVector compt_scatter_point;
  //G4double Compt_edep = 0., Photo_edep = 0.;
  // G4double KE_i = preStepPoint->GetTotalEnergy();
  // G4double KE_f = postStepPoint->GetTotalEnergy();
  
  //if (step->GetPreStepPoint()->GetProcessDefinedStep()){
  G4cout << "step pos: (" << phopos_post[0] << ", " << phopos_post[1] << ", " << phopos_post[2] << ") " << G4endl;
  G4cout << "Proc Name: " <<  proc << G4endl;
  G4cout << "step edep: " << edep << G4endl;

  G4cout << "energy at pre: " << KE_i << G4endl;

  G4cout << "energy at post: " << KE_f << G4endl;

  G4cout << "difference between pre and post: " << KE_i-KE_f << G4endl;
  //}

  if (proc == "compt" || proc == "Rayl") fEventAction-> Compt_edep.push_back(KE_i - KE_f);
  if (proc == "phot") fEventAction-> Photo_edep.push_back(KE_i - KE_f);


  //G4cout << "compton scattering postion " << compt_scatter_point << G4endl;
   
  if (proc == "eIoni")
    { 
      G4cout << "ionisation pre coord: (" << phopos_pre[0] << ", " << phopos_pre[1] << ", " << phopos_pre[2] << ") " << G4endl;
      G4cout << "ionisation post coord: (" << phopos_post[0] << ", " << phopos_post[1] << ", " << phopos_post[2] << ") " << G4endl;

      // if (phopos_pre == compt_scatter_point) {
      // 	G4cout << "energy deposited in compton: " << edep << G4endl;
      // }

      //if (step->GetPreStepPoint()->GetProcessDefinedStep())
      //G4cout << "prestep process: " << preStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;      
    }
  G4cout << G4endl;

  //G4cout << "Energy deposited via compton: " << Compt_edep << G4endl;
  //G4cout << "Energy deposited via photoelectric: " << Photo_edep << G4endl;

  G4cout << G4endl;
    
  G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  
  // if(volume != fScoringVolume)
  // return; 
  
  const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
  G4ThreeVector position = detectorConstruction->sensDet->posDetector;
  //G4ThreeVector position = detectorConstruction->sensDet->posPhoton;
  //G4double Energy_dep = detectorConstruction->sensDet->edep;

  G4double Xcord = position[0];
  G4double Ycord = position[1];
  G4double Zcord = position[2];

  
  //if (edep ==0) return;
  fEventAction->AddEdep(edep);
  //fEventAction->AddEdep(Energy_dep);
  
  fEventAction->AddX(Xcord);
  fEventAction->StorePos(Xcord, Ycord, Zcord, edep);
  
}
