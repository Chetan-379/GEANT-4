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

  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    
  //=================Studying photon interaction with material===========
  G4double edep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edep);
      
  G4double KE_i = preStepPoint->GetKineticEnergy();
  G4double KE_f = postStepPoint->GetKineticEnergy();
 
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


  //=================Analysing Scintillation photon======================
  G4String processName;
  if (track->GetCreatorProcess()) {
    processName = track->GetCreatorProcess()->GetProcessName();
  }
}
   




// G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
// const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
// G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();


