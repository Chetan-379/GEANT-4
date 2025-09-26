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
  
  // if (VolName == "physDetector" && postStepPoint->GetStepStatus() == fGeomBoundary) {
  //   track->SetTrackStatus(fStopAndKill);
  // }

  

  //G4String CreatorProc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    


  //=================Analysing Scintillation photon======================
  G4String CreatorProc;
  G4ThreeVector CreatPos = track->GetVertexPosition();
  G4ThreeVector OpPho_Pos = track->GetPosition();
  G4double OpPho_lmbda = 0, OpPho_Energy = 0;
  
  if (track->GetCreatorProcess()) {
    CreatorProc = track->GetCreatorProcess()->GetProcessName();

    if(particleName == "opticalphoton" && CreatorProc == "Scintillation"){
      G4AnalysisManager *man = G4AnalysisManager::Instance();
      OpPho_Energy = track->GetKineticEnergy();
      OpPho_lmbda = (1239.8e-6)/OpPho_Energy;   //unit: nm
      
      if(preStepPoint->GetPosition() == CreatPos){
	fEventAction->nOpPhotons++;	
	
	man->FillNtupleDColumn(0, 0, OpPho_lmbda);
	man->FillNtupleDColumn(0, 1, CreatPos[2]);
	man->AddNtupleRow(0);
      }
      
      if(OpPho_Pos[2] == 10*mm){
	fEventAction->nOpPhotons_end++;	
   
	man->FillNtupleDColumn(1, 0, OpPho_lmbda);
	man->FillNtupleDColumn(1, 1, CreatPos[2]);
	man->FillNtupleDColumn(1, 2, OpPho_Pos[2]);
	man->AddNtupleRow(1);

	track->SetTrackStatus(fStopAndKill);
      }	
    }
  }
}




   




// G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
// const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
// G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();


