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

  const G4TouchableHandle& preStepTouch = preStepPoint->GetTouchableHandle();
  const G4VPhysicalVolume* Volume = preStepTouch->GetVolume();
  const G4String& VolName = Volume->GetName();
  
  const G4VPhysicalVolume* NxtVolume = preStepTouch->GetVolume();
  const G4String& NxtVolName = Volume->GetName();

  G4Track *track = step->GetTrack();
  G4ParticleDefinition* particleDef = track->GetDefinition();
  G4String particleName = particleDef->GetParticleName();
  G4ThreeVector TrkPos = track->GetPosition();

  if (VolName == "physDetector" && postStepPoint->GetStepStatus() == fGeomBoundary) {
    track->SetTrackStatus(fStopAndKill);
  }

  G4String processName;
  if (track->GetCreatorProcess()) {
    processName = track->GetCreatorProcess()->GetProcessName();
  }

  
  G4double OptPho_Energy;
  //if (particleName == "opticalphoton" && processName == "Cerenkov") track->SetTrackStatus(fStopAndKill);
  
  if (particleName == "opticalphoton" && processName != "Cerenkov")
    {      
      fEventAction -> nOptPho++;
      fEventAction -> OptPho_PosX.push_back(TrkPos[0]);
      fEventAction -> OptPho_PosY.push_back(TrkPos[1]);
      fEventAction -> OptPho_PosZ.push_back(TrkPos[2]);
       
      if (abs(TrkPos[0]) <= 100 && abs(TrkPos[1]) <= 100 && abs(TrkPos[2]) == 150) fEventAction->TrkOnDet++;
      //G4cout << "optical photon generated from: " << processName << G4endl;
      OptPho_Energy = track->GetKineticEnergy();
	fEventAction ->OptPho_Energy.push_back(OptPho_Energy);
      
      //G4cout << "Energy of the scintillation photon is: " << OptPhoEnergy << G4endl;

      
    }
  
      // if (TrkPos[2] < -100) G4cout << "Process followed: " << proc << G4endl;
     
  
  
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

  //============
  // G4cout << "step pos: " << phopos_post << G4endl;
  // G4cout << "Proc Name: " <<  proc << G4endl;
  // G4cout << "step edep: " << edep << G4endl;

  // G4cout << "energy at pre: " << KE_i << G4endl;

  // G4cout << "energy at post: " << KE_f << G4endl;

  // G4cout << "difference between pre and post: " << KE_i-KE_f << G4endl;

  // G4cout << "parentID: " << step-> GetTrack()->GetParentID() << G4endl;
  //========================
  //}

  if (step->GetTrack()->GetParentID() == 0){
    if (proc == "compt" || proc == "Rayl") {
      fEventAction-> Compt_edep.push_back(KE_i - KE_f);
    }
    
    if (proc == "phot") {
      
      fEventAction-> Photo_edep.push_back(KE_i - KE_f);
      // if (fEventAction->Photo_edep.size() == 0 && fEventAction->Compt_edep.size() == 0) fEventAction->Compt_Before_Photo = false;
      // if (fEventAction->Photo_edep.size() == 0 && fEventAction->Compt_edep.size() > 0) fEventAction->Compt_Before_Photo = true;
      // if (fEventAction->Photo_edep.size() == 0) fEventAction->nCompton_Before_Photo = Compt_edep.size();
      
      //break;
    }
    fEventAction->Xarray.push_back(phopos_post[0]);
    fEventAction->Yarray.push_back(phopos_post[1]);
    fEventAction->Zarray.push_back(phopos_post[2]);
  }
    
  
  //G4cout << "compton scattering postion " << compt_scatter_point << G4endl;
   
  if (proc == "eIoni")
    { //====================
      // G4cout << "ionisation pre coord: " << phopos_pre << G4endl;
      // G4cout << "ionisation post coord: " << phopos_post << G4endl;
      //===================

      
      //if (step->GetPreStepPoint()->GetProcessDefinedStep())
      //G4cout << "prestep process: " << preStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;      
    }
  //G4cout << G4endl;

  // G4cout << "Energy deposited via compton: " << Compt_edep << G4endl;
  // G4cout << "Energy deposited via photoelectric: " << Photo_edep << G4endl;

  // G4cout << G4endl;


    
  G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
  // if(volume != fScoringVolume)
  // return; 
  
  const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
  //G4ThreeVector position = detectorConstruction->sensDet->posDetector;
  G4ThreeVector position = phopos_post;
  //G4ThreeVector position = detectorConstruction->sensDet->posPhoton;
  //G4double Energy_dep = detectorConstruction->sensDet->edep;

  G4double Xcord = position[0];
  G4double Ycord = position[1];
  G4double Zcord = position[2];
  
  fEventAction->AddEdep(edep);
  fEventAction->AddX(Xcord);
  //fEventAction->StorePos(Xcord, Ycord, Zcord, edep);
  
}
