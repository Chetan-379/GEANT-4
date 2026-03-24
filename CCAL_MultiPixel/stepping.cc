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

  

  G4double edep = step->GetTotalEnergyDeposit();


  //=================Analysing Scintillation photon======================
  G4TouchableHistory *postTouch = (G4TouchableHistory *)(postStepPoint->GetTouchable());
  G4TouchableHistory *preTouch = (G4TouchableHistory *)(preStepPoint->GetTouchable());

  G4String VolName = preTouch->GetVolume()->GetName();
  //G4String VolName = track->GetTouchable()->GetVolume()->GetName();

  G4String creatProc;
  if (track->GetCreatorProcess())    creatProc = track->GetCreatorProcess()->GetProcessName();
 
  G4float Op_lmbda = (1239.8e-6) / (postStepPoint->GetKineticEnergy());

  if (VolName == "GAGG_PV" || VolName == "LYSO_PV")
    {
      G4int Det_cpyNo = preStepPoint->GetTouchable()->GetVolume()->GetCopyNo();
      //G4int Det_cpyNo = track->GetTouchable()->GetVolume()->GetCopyNo();
      G4int clr_Idx = abs(Det_cpyNo/100);
      G4int row_Idx = abs(G4int ((Det_cpyNo%100)/10))-1;
      G4int col_Idx = abs((Det_cpyNo%100) % 10)-1;

      G4String proc;
      auto Process = step->GetPostStepPoint()->GetProcessDefinedStep();
      if(Process != nullptr) proc = Process->GetProcessName();
  
      if (particleName == "opticalphoton" && creatProc == "Scintillation")
	{
	  bool isGAGG = false, isLYSO = false;       
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
	  
	  if (nullptr == fBoundary) {
	    G4ProcessManager* pm = part->GetProcessManager();
	    G4int nprocesses = pm->GetProcessListLength();
	    G4ProcessVector* pv = pm->GetProcessList();
	    for (G4int i = 0; i < nprocesses; ++i) {
	      if (nullptr != (*pv)[i] && (*pv)[i]->GetProcessName() == "OpBoundary") {
		fBoundary = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
		break;  // find the boundary process only once
	      }
	    }
	  }
      
	  if (nullptr != fBoundary) boundaryStatus = fBoundary->GetStatus();    
	  if (boundaryStatus == Detection) {   //detection only possible at SiPM surface as other have zero efficiency
	    fRunAction->h_Op_QE_lmbda->Fill(Op_lmbda);
	    fRunAction->h_Op_QE_energy->Fill(postStepPoint->GetKineticEnergy()*1e6);
	    fRunAction->h_Op_QE_Z->Fill(TrkPos[2]);

	    if(isGAGG) fRunAction->nOpGAGG_QE++;
	    if(isLYSO) fRunAction->nOpLYSO_QE++;
	  }
	}
      
      fRunAction->Edep_truth += edep;      
      if (VolName == "GAGG_PV") fRunAction->Edep_G_truth +=edep;	            
      if (VolName == "LYSO_PV") fRunAction->Edep_L_truth +=edep;	      
      
      if (edep > 0.00001) {
	bool updateEntry = false, addEntry =false;
	G4int updt_Idx = -999;
	G4int nCompt =0, nPhoto =0;
	G4double scat_theta = -1000;
	G4ThreeVector vin = preStepPoint->GetMomentumDirection();
	G4ThreeVector vout = postStepPoint->GetMomentumDirection();

	if(parentId == 0 && proc == "compt") {
	  nCompt++;
	  //scat_theta = acos(vin.dot(vout));

	  auto secondaries = step->GetSecondaryInCurrentStep();
	  if (secondaries->size()>0){
	    auto scat_e_E = (*secondaries)[0]->GetKineticEnergy();	 	 
	    
	    scat_theta = acos(1-(0.511*((1/(0.511-scat_e_E)) - (1/0.511))));
	    //G4cout << "scat theta: " << scat_theta << G4endl;
	  }
	}
	
	if(parentId == 0 && proc == "phot") nPhoto++;
	
	for(int i =0; i< fRunAction->Det_row_Idx.size(); i++){
	  G4int check_DetId = 100*(fRunAction->Det_clr_Idx[i]) + 10*(fRunAction->Det_row_Idx[i]+1) + fRunAction->Det_col_Idx[i] +1;	  
	  
	  if(Det_cpyNo == check_DetId) {
	    updateEntry =true;
	    addEntry =false;
	    updt_Idx = i;	    
	    break;
	  }
	  
	  else addEntry =true;       
	}
	
	if (addEntry){
	  fRunAction->Det_row_Idx.push_back(row_Idx);
	  fRunAction->Det_col_Idx.push_back(col_Idx);
	  fRunAction->Det_clr_Idx.push_back(clr_Idx);
	  fRunAction->Edep_truth_vec.push_back(edep);
	  fRunAction->Det_nCompt.push_back(nCompt);
	  fRunAction->Det_nPhoto.push_back(nPhoto);
	  fRunAction->Det_theta.push_back(scat_theta);
	}
	
	if (updateEntry) {
	  fRunAction->Edep_truth_vec[updt_Idx] += edep;
	  fRunAction->Det_nCompt[updt_Idx] += nCompt;
	  fRunAction->Det_nPhoto[updt_Idx] += nPhoto;
	}

	if (fEventAction->storeHit) {
	  fRunAction->Det_row_Idx.push_back(row_Idx);
	  fRunAction->Det_col_Idx.push_back(col_Idx);
	  fRunAction->Det_clr_Idx.push_back(clr_Idx);
	  fRunAction->Edep_truth_vec.push_back(edep);	  
	  fRunAction->Det_nCompt.push_back(nCompt);
	  fRunAction->Det_nPhoto.push_back(nPhoto);
	  fRunAction->Det_theta.push_back(scat_theta);

	  fEventAction->storeHit = false; //only for storing the 1st interaction of an event
	}  
      }      
    }  
}

