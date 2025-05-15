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
  G4ThreeVector TrkPos = track->GetPosition();

  G4ParticleDefinition* particleDef = track->GetDefinition();
  G4String particleName = particleDef->GetParticleName();
  
  if (VolName == "physDetector" && postStepPoint->GetStepStatus() == fGeomBoundary) {
    track->SetTrackStatus(fStopAndKill);
  }

  G4String processName;
  if (track->GetCreatorProcess()) {
    processName = track->GetCreatorProcess()->GetProcessName();
  }
  
  G4double OptPho_Energy;
  G4double OptPhot_time;

  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  
  if (particleName == "opticalphoton" && processName != "Cerenkov")
    {
      fEventAction -> nOptPho++;
      fEventAction -> OptPho_PosX.push_back(TrkPos[0]);
      fEventAction -> OptPho_PosY.push_back(TrkPos[1]);
      fEventAction -> OptPho_PosZ.push_back(TrkPos[2]);
      
      OptPho_Energy = track->GetKineticEnergy();
      fEventAction ->OptPho_Energy.push_back(OptPho_Energy);

      OptPhot_time = track->GetGlobalTime();
      fEventAction ->OptPho_time.push_back(OptPhot_time);
      if (abs(TrkPos[0]) <= 100 && abs(TrkPos[1]) <= 100 && abs(TrkPos[2]) == 150) fEventAction->TrkOnDet++;
    }
  

  G4double edep = step->GetTotalEnergyDeposit();
  
    
  G4double KE_i = preStepPoint->GetKineticEnergy();
  G4double KE_f = postStepPoint->GetKineticEnergy();

  G4ThreeVector compt_scatter_point;

 
  if (step->GetTrack()->GetParentID() == 0){
    if (proc == "compt" || proc == "Rayl") {
      fEventAction-> Compt_edep.push_back(KE_i - KE_f);
    }
    
    if (proc == "phot") {
      
      fEventAction-> Photo_edep.push_back(KE_i - KE_f);
    }
    fEventAction->Xarray.push_back(phopos_post[0]);
    fEventAction->Yarray.push_back(phopos_post[1]);
    fEventAction->Zarray.push_back(phopos_post[2]);
  }
          
  G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
  const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
  G4ThreeVector position = phopos_post;
 
  G4double Xcord = position[0];
  G4double Ycord = position[1];
  G4double Zcord = position[2];
  
  fEventAction->AddEdep(edep);
  fEventAction->AddX(Xcord); 
}

