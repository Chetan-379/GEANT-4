#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction, MyRunAction* runAction)
{
  fEventAction = eventAction;
  fRunAction = runAction;
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
  
  G4Track *track = step->GetTrack();

  G4int parentId = track->GetParentID();
  G4ThreeVector TrkPos = track->GetPosition();

  G4int TrkId = track-> GetTrackID();

  G4ParticleDefinition* particleDef = track->GetDefinition();
  G4String particleName = particleDef->GetParticleName();
  G4ThreeVector creatPos = track->GetVertexPosition();
  
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
 
  //=================Studying photon interaction with material===========
  G4double edep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edep);
  
  G4double KE_i = preStepPoint->GetKineticEnergy();
  G4double KE_f = postStepPoint->GetKineticEnergy();
  
  if (step->GetTrack()->GetParentID() == 0 && TrkId ==1 && TrkPos[2]>0){
    if (proc == "compt") {
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
  
  G4TouchableHistory* postTouch = (G4TouchableHistory*) (postStepPoint->GetTouchable());
  G4TouchableHistory* preTouch = (G4TouchableHistory*) (preStepPoint->GetTouchable());
  
  G4String VolName= preTouch->GetVolume()->GetName();

  auto history = preTouch->GetHistory();
  G4int depth = history->GetDepth(); 
  
  G4String creatProc;
  if (track->GetCreatorProcess()) creatProc = track->GetCreatorProcess()->GetProcessName();
  
  G4int DetId = 999;  
  G4ThreeVector prePos = preStepPoint->GetPosition();

  G4ThreeVector DetPos;
 
  if (VolName == "Scint_PV"){
    DetId = (10*(preTouch->GetVolume(1)->GetCopyNo()+1)+preTouch->GetVolume(2)->GetCopyNo()+1) * preTouch->GetVolume(3)->GetCopyNo();   
    
    G4int row_Idx, col_Idx, module_Idx;
    
    row_Idx = preTouch->GetVolume(1)->GetCopyNo();
    col_Idx = preTouch->GetVolume(2)->GetCopyNo();
    module_Idx = preTouch->GetVolume(3)->GetCopyNo();
    
    DetPos =  preTouch->GetTranslation();

    // if (edep>0) G4cout << module_Idx << "*[" << row_Idx << "][" << col_Idx << "]" << G4endl;

    if(particleName == "opticalphoton" && creatProc == "Scintillation"){
      if (prePos == creatPos) {
	if (module_Idx == 1) {
	  fEventAction->matrix1_gen[row_Idx][col_Idx]++;
	  fRunAction->Module1[row_Idx][col_Idx][nGenOp]++;	  
	}
      
	if (module_Idx == -1){
	  fEventAction->matrix2_gen[row_Idx][col_Idx]++;
	  fRunAction->Module1[row_Idx][col_Idx][nGenOp]++;
	}
      }
    
      if (TrkPos[2] == 35*mm || TrkPos[2] == -35*mm){
	if (module_Idx == 1) {
	  fEventAction->matrix1[row_Idx][col_Idx]++;
	
	  fRunAction->Module1[row_Idx][col_Idx][nOpPho]++;
	  fRunAction->Module1[row_Idx][col_Idx][DetPosX] = DetPos[0];
	  fRunAction->Module1[row_Idx][col_Idx][DetPosY] = DetPos[1];
	}
      
	if (module_Idx == -1){
	  fEventAction->matrix2[row_Idx][col_Idx]++;
	
	  fRunAction->Module2[row_Idx][col_Idx][nOpPho]++;
	  fRunAction->Module2[row_Idx][col_Idx][DetPosX] = DetPos[0];
	  fRunAction->Module2[row_Idx][col_Idx][DetPosY] = DetPos[1];	  
	}
      
	//G4cout << "DetPos: " << DetPos << G4endl;
      
	track->SetTrackStatus(fStopAndKill);   //killing the Op pho at rear surface
      }
    }

    else if (particleName != "opticalphoton"){
      if(module_Idx == 1) fRunAction->Module1[row_Idx][col_Idx][Edep] += edep;
      if(module_Idx == -1) fRunAction->Module2[row_Idx][col_Idx][Edep] += edep;      
    }      
  }
}
           
  
   

   




// G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  
// const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
// G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();


// for (int i = 0; i < depth; ++i) {
// 	G4cout << "[" << i << "] "
// 	       << history->GetVolume(i)->GetName()
// 	       << " copy=" << history->GetReplicaNo(i)
// 	       << G4endl;
// }
