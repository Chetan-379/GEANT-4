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

  //const G4TouchableHandle& preStepTouch = preStepPoint->GetTouchableHandle();
  const G4VPhysicalVolume* PreVolume = preStepPoint->GetTouchableHandle()->GetVolume();
  const G4VPhysicalVolume* PostVolume = postStepPoint->GetTouchableHandle()->GetVolume();
  
  const G4String& PreVolName = PreVolume->GetName();
  const G4String& PostVolName = PostVolume->GetName();
  
  // const G4VPhysicalVolume* NxtVolume = preStepTouch->GetVolume();
  // const G4String& NxtVolName = Volume->GetName();

  G4Track *track = step->GetTrack();
  G4ThreeVector TrkPos = track->GetPosition();

  G4ParticleDefinition* particleDef = track->GetDefinition();
  G4String particleName = particleDef->GetParticleName();
  
  // if (PreVolName == "physDetector" && postStepPoint->GetStepStatus() == fGeomBoundary) {
  //   //track->SetTrackStatus(fStopAndKill);
  // }

  G4String processName;
  if (track->GetCreatorProcess()) {
    processName = track->GetCreatorProcess()->GetProcessName();
  }
  
  G4double OptPho_Energy;
  G4double OptPhot_time;

  // if (particleName == "opticalphoton" && processName == "Cerenkov") {
  //   track->SetTrackStatus(fStopAndKill);
  //   //G4cout << "creator process is: " << processName << G4endl;
  // }

   // if (particleName == "opticalphoton" && PreVolName == "physDetector" && PostVolName == "physEncVol") {
   //   G4cout << "optical photon moved out" << G4endl;

     //track->SetTrackStatus(fStopAndKill);
  //   }     
  // }

  // if (particleName == "opticalphoton" && track->GetVolume()->GetName() == "physEncVol") {G4cout << "optical photon moved out" << G4endl;
  //   track->SetTrackStatus(fStopAndKill);}
   
  //G4int nBoundProc = 0, nRefraction =0, nReflection = 0, nTIR =0;
  if (particleName == "opticalphoton" && processName != "Cerenkov")
    {  
      fEventAction -> nOptPho++;
      fEventAction -> OptPho_PosX.push_back(TrkPos[0]);
      fEventAction -> OptPho_PosY.push_back(TrkPos[1]);
      fEventAction -> OptPho_PosZ.push_back(TrkPos[2]);
       
      if (abs(TrkPos[0]) <= 100 && abs(TrkPos[1]) <= 100 && abs(TrkPos[2]) == 150) fEventAction->TrkOnDet++;

      OptPho_Energy = track->GetKineticEnergy();
      fEventAction ->OptPho_Energy.push_back(OptPho_Energy);

      OptPhot_time = track->GetGlobalTime();
      fEventAction ->OptPho_time.push_back(OptPhot_time);
      
      if (postStepPoint->GetStepStatus() == fGeomBoundary) {
      	G4OpBoundaryProcessStatus theStatus = Undefined;
      	G4ProcessManager* opManager = G4OpticalPhoton::OpticalPhotonDefinition()->GetProcessManager();
      	G4int n_proc = opManager->GetPostStepProcessVector(typeDoIt)->entries();
      	G4ProcessVector* postStepDoItVector = opManager->GetPostStepProcessVector(typeDoIt);

      	//G4cout << "total no. of boundary processes occured: " << n_proc << G4endl;
      	for (G4int i = 0; i < n_proc; ++i) {
      	  G4VProcess* currentProcess = (*postStepDoItVector)[i];
      	  auto opProc = dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
      	  if (opProc) theStatus = opProc->GetStatus();
      	  //G4cout << "\n\nOptical process at boundary is: " << theStatus << G4endl;

      	  if (theStatus != Undefined && theStatus != NotAtBoundary && theStatus != StepTooSmall) {
      	  switch (theStatus) {
	      
      	  case FresnelReflection:
      	    G4cout << "Fresnel Reflection" << G4endl;
      	    fEventAction ->n_Reflection++;
      	    fEventAction ->n_BoundProc++;
      	    break;

      	  case FresnelRefraction:
      	    //G4cout << "Fresnel Refraction" << G4endl;
      	    fEventAction ->n_Refraction++;
      	    fEventAction ->n_BoundProc++;
	    fEventAction ->nOpPhoLeft++;
      	    break;

      	  case TotalInternalReflection:
      	    G4cout << "Total Internal Reflection" << G4endl;
      	    fEventAction ->n_TIR++;
      	    fEventAction ->n_BoundProc++; 
      	    break;

      	  case LobeReflection:
      	    G4cout << "Lobe Reflection" << G4endl;
      	    break;

      	  case BackScattering:
      	    //G4cout << "Back Scattering" << G4endl;
      	    break;

      	  case Absorption:
      	    //G4cout << "Absorption" << G4endl;
      	    break;

      	  case Detection:
      	    //G4cout << "Detection" << G4endl;
      	    break;
     
      	  case Transmission:
      	    //G4cout << "Transmission" << G4endl;
      	    break;

      	  default:
      	    //G4cout << "Other Optical Boundary Status" << G4endl;
      	    break;
      	  }

      	  // fEventAction->n_Refraction = nRefraction;      
      	  // fEventAction->n_Reflection = nReflection;
      	  // fEventAction->n_TIR = nTIR;
	  
      	  // //if (nBoundProc >1){
      	  //   G4cout << "no. of boundary process: " << nBoundProc << G4endl;    
      	  //   G4cout << "no. of refraction: " << nRefraction << G4endl;
      	  //   G4cout << "no. of reflection: " << nReflection << G4endl;
      	  //   G4cout << "no. of TIR: " << nTIR << G4endl;}
      	  // // }
	
      	// fEventAction->AddBoundary();
      	  }
      	}
      }
    }

      
 
      
  G4double edep = step->GetTotalEnergyDeposit();
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    
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

  
}

