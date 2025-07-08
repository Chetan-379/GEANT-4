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

  //getting the position of compton scattering point
  //G4String proc = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4ThreeVector compt_scat_point;
  G4int ParentID = step->GetTrack()->GetParentID();
  
  // G4int nCompt = 0;
  // if(proc == "compt" && ParentID == 0) nCompt++;
  // G4cout << "\n nCompt: " << nCompt << G4endl;
    
  
  // Get hit accounting data for this cell
  auto hit = new CellHit();

  // Add values
  hit->Add(edep, stepLength);
 
  //set position, time etc
  G4double currentTime = step->GetPreStepPoint()->GetGlobalTime();
  G4double TrkLen = step->GetTrack()->GetTrackLength();
  G4int DetId = touchable->GetVolume()->GetCopyNo();
  G4double scat_theta = -1;
  G4ThreeVector vin = step->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vout = step->GetPostStepPoint()->GetMomentumDirection();
  G4double eta = 10000;
  
  // G4ThreeVector vin = step->GetPreStepPoint()->GetMomentum();
  // G4ThreeVector vout = step->GetPostStepPoint()->GetMomentum();

  // G4cout << "vin mag" << vin.mag() << G4endl;

  if ((proc == "compt" || proc == "phot" || proc == "Rayl") && ParentID == 0) {
    compt_scat_point = step->GetPostStepPoint()->GetPosition();

    G4cout << "\nCompton Scat point: " << compt_scat_point << G4endl;
    G4cout << "pre step momDir: " << vin << G4endl;
    G4cout << "post step momDir: " << vout << G4endl;

    if (proc != "phot") {scat_theta = acos(vin.dot(vout));
      auto prePol = step->GetPreStepPoint()->GetPolarization();
      auto postPol = step->GetPostStepPoint()->GetPolarization();

    G4cout << "polarisation before scattering: " << prePol << G4endl;
    G4cout << "polarisation after scattering: " << postPol << "\n\n";

    auto scatPlane = vin.cross(vout);
    auto polPlane = vin.cross(prePol);
    //autoq polPlane = vin;

    
    //auto prfY = G4PolarizationHelper::GetParticleFrameY(vin);
    auto prfY = G4PolarizationHelper::GetParticleFrameY(vin);
    auto prfX = G4PolarizationHelper::GetParticleFrameX(vin);

    //G4cout << "dot of vin with pol: " << vin.dot(prf) << "\n\n";

    //G4PolarizationHelper::TestPolarizationTransformations();
    //auto post_prf = G

    G4cout << "\nvin: " << vin << "\n";
    G4cout << "PRF_Xaxis:" << prfX << "\n";
    G4cout << "PRF_Yaxis:" << prfY << "\n";
    
    //eta = acos((scatPlane.dot(polPlane))/(scatPlane.mag()*polPlane.mag()));
    auto vinPol_Angle = vin.dot(prfY);

    // if (vinPol_Angle > 0.000001) G4cout << "\n\n\n=============================ALERT!!==================" << G4endl;

    G4cout << "\nvin.prfY: " << vinPol_Angle << G4endl;
    eta = acos((scatPlane.dot(prfY))/(scatPlane.mag()*prfY.mag()));

    //G4cout << "\n\neta is: " << eta << G4endl;
    }
    
    hit->SetPosition(compt_scat_point);
    hit->SetTime(currentTime);
    hit->SetDetectorID(DetId);
    hit->SetTrackLength(TrkLen);
    hit->SetProcName(proc);
    hit->SetProcId(procId);
    hit->SetScatAngle(scat_theta);
    hit->SetEta(eta);

    fHitsCollection->insert(hit);

    //getting the polarisation of the photon after different compton scattering:

    
    
  }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    auto nofHits = fHitsCollection->entries();
    
    G4cout << "\n" << "-------->Hit summary of the SD "  << G4endl;
    for (std::size_t i = 0; i < nofHits; ++i)
      (*fHitsCollection)[i]->Print();
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
