#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction, MyRunAction *runAction)
{
  fEventAction = eventAction;
  fRunAction = runAction;
}

MySteppingAction::~MySteppingAction()
{
}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{
  G4StepPoint *preStepPoint = step->GetPreStepPoint();
  G4StepPoint *postStepPoint = step->GetPostStepPoint();

  G4ThreeVector phopos_pre = preStepPoint->GetPosition();
  G4ThreeVector phopos_post = postStepPoint->GetPosition();

  const G4TouchableHandle &preStepTouch = preStepPoint->GetTouchableHandle();

  G4Track *track = step->GetTrack();

  G4int parentId = track->GetParentID();
  G4ThreeVector TrkPos = track->GetPosition();

  G4int TrkId = track->GetTrackID();

  G4ParticleDefinition *particleDef = track->GetDefinition();
  G4String particleName = particleDef->GetParticleName();
  G4ThreeVector creatPos = track->GetVertexPosition();

  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  //=================Studying photon interaction with material===========
  G4double edep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edep);

  G4double KE_i = preStepPoint->GetKineticEnergy();
  G4double KE_f = postStepPoint->GetKineticEnergy();

  if (step->GetTrack()->GetParentID() == 0 && TrkId == 1 && TrkPos[2] > 0)
  {
    if (proc == "compt")
    {
      fEventAction->Compt_edep.push_back(KE_i - KE_f);
    }

    if (proc == "phot")
    {
      fEventAction->Photo_edep.push_back(KE_i - KE_f);
    }

    fEventAction->Xarray.push_back(phopos_post[0]);
    fEventAction->Yarray.push_back(phopos_post[1]);
    fEventAction->Zarray.push_back(phopos_post[2]);
  }

  //=================Analysing Scintillation photon======================

  G4TouchableHistory *postTouch = (G4TouchableHistory *)(postStepPoint->GetTouchable());
  G4TouchableHistory *preTouch = (G4TouchableHistory *)(preStepPoint->GetTouchable());

  G4String VolName = preTouch->GetVolume()->GetName();

  G4String creatProc;
  if (track->GetCreatorProcess())    creatProc = track->GetCreatorProcess()->GetProcessName();
 
  G4float Op_lmbda = (1239.8e-6) / (postStepPoint->GetKineticEnergy());
  if (VolName == "GAGG_PV" || VolName == "LYSO_PV")
  {
    if (particleName == "opticalphoton" && creatProc == "Scintillation")
    {      
      if (TrkPos[2] == 35 * mm)
        {
	  G4float OpX = TrkPos[0];
	  G4float OpY = TrkPos[1];

	  G4int quad_idx = 9999;
	  if(OpX >=0 && OpY >=0) quad_idx = 0;
	  if(OpX <0 && OpY >0) quad_idx = 1;
	  if(OpX <=0 && OpY <=0) quad_idx = 2;
	  if(OpX >0 && OpY <0) quad_idx = 3;

	  fRunAction->h_zVslmbda->Fill(creatPos[2], Op_lmbda);
	  fRunAction->h_Op_lmbda->Fill(Op_lmbda);
	  fRunAction->h_Op_energy->Fill(postStepPoint->GetKineticEnergy()*1e6);
	  

	  //discrimination based on Op wavelength-------------------
	  bool isGAGG = false, isLYSO = false;       
	  if(Op_lmbda >= 480) isGAGG = true;
	  else if(Op_lmbda <= 479) isLYSO = true;
	  
	  if(isGAGG){
	    fRunAction->nOpGAGG++;
	    fRunAction->nOpG_quad[quad_idx]++;	    
	  }

	  if(isLYSO){
	    fRunAction->nOpLYSO++;
	    fRunAction->nOpL_quad[quad_idx]++;	  
	  }	  

	  //discrimination based on creation positionZ--------------
	  bool isGAGG_truth = false, isLYSO_truth = false;
	  if(creatPos[2] <= 25) isGAGG_truth =true;
	  else isLYSO_truth = true;

	  if(isGAGG_truth) fRunAction->nOpGAGG_truth++;
	  if(isLYSO_truth) fRunAction->nOpLYSO_truth++;
	  //track->SetTrackStatus(fStopAndKill); // killing the Op pho at rear surface
        }

      //counting the Op photons which gets detected after QE of SiPM
      G4OpBoundaryProcessStatus boundaryStatus = Undefined;

      auto part = track->GetDefinition();
      // find the boundary process only once
      if (nullptr == fBoundary) {
	G4ProcessManager* pm = part->GetProcessManager();
	G4int nprocesses = pm->GetProcessListLength();
	G4ProcessVector* pv = pm->GetProcessList();
	for (G4int i = 0; i < nprocesses; ++i) {
	  if (nullptr != (*pv)[i] && (*pv)[i]->GetProcessName() == "OpBoundary") {
	    fBoundary = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
	    break;
	  }
	}
      }
      
      if (nullptr != fBoundary) boundaryStatus = fBoundary->GetStatus();    
      if (boundaryStatus == Detection) {   //detection only possible at SiPM surface as other have zero efficiency
	fRunAction->h_Op_QE_lmbda->Fill(Op_lmbda);
	fRunAction->h_Op_QE_energy->Fill(postStepPoint->GetKineticEnergy()*1e6);
	fRunAction->h_Op_QE_Z->Fill(TrkPos[2]);
      }     
    }
         
    fRunAction->Edep_truth += edep;
    if (VolName == "GAGG_PV") fRunAction->Edep_G_truth +=edep;
    if (VolName == "LYSO_PV") fRunAction->Edep_L_truth +=edep;
  }  
}

