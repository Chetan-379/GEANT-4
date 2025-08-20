//taken from G4 example basic/B4/B4c

#include "CellSD.hh"

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellSD::CellSD(const G4String& name, const G4String& hitsCollectionName,
	       G4int nofCells)
  : G4VSensitiveDetector(name), fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new CellHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CellSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  //bool summary = false; 
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4int procId;

  //assigning variables to the processes
  if (proc == "Rayl") procId = 0;
  if (proc == "phot") procId = 1;
  if (proc == "compt") procId = 2;

  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
    stepLength = step->GetStepLength();
  }

  if (edep == 0. && stepLength == 0. && proc != "Rayl") return false;
  //if (stepLength == 0) return false;

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  G4ThreeVector compt_scat_point;    
  G4int ParentID = step->GetTrack()->GetParentID();
        
  // Get hit accounting data for this cell
  auto hit = new CellHit();

  // Add values
  hit->Add(edep, stepLength);
 
  //set position, time etc
  G4double currentTime = step->GetPreStepPoint()->GetGlobalTime();
  G4double TrkLen = step->GetTrack()->GetTrackLength();
  G4int DetId = touchable->GetVolume()->GetCopyNo();
  G4double scat_theta = -1000;
  G4ThreeVector vin = step->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vout = step->GetPostStepPoint()->GetMomentumDirection();
  G4double eta = 10000;
  G4ThreeVector prePol(-1000,-1000,-1000);
  G4double Ein = step->GetPreStepPoint()->GetKineticEnergy();
  G4double Eout = step->GetPostStepPoint()->GetKineticEnergy();
  
  

  if ((proc == "compt" || proc == "phot" || proc == "Rayl") && ParentID == 0) {
    compt_scat_point = step->GetPostStepPoint()->GetPosition();

    //G4cout << "\t" << proc << " Hit Pos: " << compt_scat_point << G4endl;

    if (proc != "phot") {      
      //if(step->GetPreStepPoint()->GetProcessDefinedStep() && step->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() != "Rayl"){
      //if(proc == "compt" && step->GetPreStepPoint()->GetKineticEnergy() == 0.511) scat_theta = acos(vin.dot(vout));
      //}

      scat_theta = acos(vin.dot(vout));
      
      prePol = step->GetPreStepPoint()->GetPolarization();
      auto postPol = step->GetPostStepPoint()->GetPolarization();

      
      auto scatPlane = vin.cross(vout);
      auto polPlane = vin.cross(prePol);
      auto postPolPlane = vout.cross(postPol);

      // auto prfY_before = G4PolarizationHelper::GetParticleFrameY(vin);
      // auto prfX_before = G4PolarizationHelper::GetParticleFrameX(vin);

      auto prfY_after = G4PolarizationHelper::GetParticleFrameY(vout);
      auto prfX_after = G4PolarizationHelper::GetParticleFrameX(vout);

      // float pol_mom_Angle_bef = acos(vin.dot(prfX_before)) * 180/CLHEP::pi;
      // float pol_mom_Angle_aft = acos(vout.dot(prfX_after)) * 180/CLHEP::pi;

      float pol_mom_Angle_bef = acos(vin.dot(prePol)) * 180/CLHEP::pi;
      float pol_mom_Angle_aft = acos(vout.dot(postPol)) * 180/CLHEP::pi;

      //================ working eta===============
      //eta = std::atan2(vin.dot(prePol.cross(vout)), prePol.dot(vout));
      //===========================================

      eta = vout.azimAngle(prePol,vin);
      //eta = vout.azimAngle(polPlane,vin);

      //eta = acos((scatPlane.dot(polPlane))/(scatPlane.mag()*polPlane.mag()));
      //eta = acos((scatPlane.dot(prfX_before))/(scatPlane.mag()*prfX_before.mag()));
      //eta = std::atan2(vin.dot(prfX_before.cross(vout)), prfX_before.dot(vout));
      
      if(summary){  
      G4cout << "polarization vector before scattering: " << prePol << G4endl;
      G4cout << "polarization vector after scattering: " << postPol << "\n\n";
      
      G4cout << "vin Direction: " << vin << "\n";
      G4cout << "polarisation (Elec field Direction) before scattering: " << prePol << "\n";
      G4cout << "polPlane normal before scattering: " << polPlane << "\n\n";
      
      G4cout << "vout Direction: " << vout << "\n";
      G4cout << "polarisation (Elec field Direction) after scattering:" << postPol << "\n";
      G4cout << "polPlane normal after scattering:" << postPolPlane << "\n\n";
            
      G4cout << "Angle between polDir & ptcl momDir before scattering: " << pol_mom_Angle_bef << " deg" << G4endl;
      G4cout << "Angle between polDir & ptcl momDir after scattering: " << pol_mom_Angle_aft << " deg" << "\n\n";
    
      G4cout << "Angle between polPlane & ScatPlane: " << eta*180/CLHEP::pi  << " deg" << "\n\n";

      G4cout << "photon Energy before Scattering: " << Ein << "\n";
      G4cout << "photon Energy after Scattering: " << Eout << "\n";

      
      G4cout << "\n*****" << "\n\n";
      }

      // if (proc == "compt" && scat_theta*180/CLHEP::pi > 89 && scat_theta*180/CLHEP::pi < 91){
      // G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();\

      // G4cout << "\n\n***************Event with 90 deg compt: " << eventID << G4endl;}
    }
    
    hit->SetPosition(compt_scat_point);
    hit->SetTime(currentTime);
    hit->SetDetectorID(DetId);
    hit->SetTrackLength(TrkLen);
    hit->SetProcName(proc);
    hit->SetProcId(procId);
    hit->SetScatAngle(scat_theta);
    hit->SetEta(eta);
    hit->SetPolarisation(prePol);
    hit->SetEin(Ein);
    hit->SetEout(Eout);

    fHitsCollection->insert(hit);
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1 && summary) {
    auto nofHits = fHitsCollection->entries();
    
    G4cout << "\n" << "-------->Hit summary of the SD "  << G4endl;
    for (std::size_t i = 0; i < nofHits; ++i)
      (*fHitsCollection)[i]->Print();
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
