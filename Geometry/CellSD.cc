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
  G4double eta = 10000, eta1 = 10000;
  G4ThreeVector prePol(-1000,-1000,-1000);
  
  // G4ThreeVector vin = step->GetPreStepPoint()->GetMomentum();
  // G4ThreeVector vout = step->GetPostStepPoint()->GetMomentum();

  // G4cout << "vin mag" << vin.mag() << G4endl;

  if ((proc == "compt" || proc == "phot" || proc == "Rayl") && ParentID == 0) {
    compt_scat_point = step->GetPostStepPoint()->GetPosition();

    G4cout << "\t" << proc << " Hit Pos: " << compt_scat_point << G4endl;
    // G4cout << "pre step momDir: " << vin << G4endl;
    // G4cout << "post step momDir: " << vout << "\n\n";

    if (proc != "phot") {scat_theta = acos(vin.dot(vout));
      prePol = step->GetPreStepPoint()->GetPolarization();
      auto postPol = step->GetPostStepPoint()->GetPolarization();

    G4cout << "stokes vector before scattering: " << prePol << G4endl;
    G4cout << "stokes vector after scattering: " << postPol << "\n\n";

    auto scatPlane = vin.cross(vout);
    auto polPlane = vin.cross(prePol);
    //autoq polPlane = vin;

    //auto prfY = G4PolarizationHelper::GetParticleFrameY(vin);
    auto prfY_before = G4PolarizationHelper::GetParticleFrameY(vin);
    auto prfX_before = G4PolarizationHelper::GetParticleFrameX(vin);

    auto prfY_after = G4PolarizationHelper::GetParticleFrameY(vout);
    auto prfX_after = G4PolarizationHelper::GetParticleFrameX(vout);


    float pol_mom_Angle_bef = acos(vin.dot(prfX_before)) * 180/CLHEP::pi;
    float pol_mom_Angle_aft = acos(vout.dot(prfX_after)) * 180/CLHEP::pi;

    //G4PolarizationHelper::TestPolarizationTransformations();
    //auto post_prf = G

    G4cout << "vin Direction: " << vin << "\n";
    G4cout << "polarisation Direction (PRF_X) before scattering: " << prfX_before << "\n";
    G4cout << "polPlane normal (PRF_Y) before scattering: " << prfY_before << "\n\n";

    G4cout << "vout Direction: " << vout << "\n";
    G4cout << "polarisation Direction (PRF_X) after scattering:" << prfX_after << "\n";
    G4cout << "polPlane normal (PRF_Y) after scattering:" << prfY_after << "\n\n";
    
    //eta = acos((scatPlane.dot(polPlane))/(scatPlane.mag()*polPlane.mag()));
    //auto vinPol_Angle = vin.dot(prfY_before);

    // if (vinPol_Angle > 0.000001) G4cout << "\n\n\n=============================ALERT!!==================" << G4endl;

    //G4cout << "\nvin.prfY: " << vinPol_Angle << G4endl;

    G4cout << "Angle between polDir & ptcl momDir before scattering: " << pol_mom_Angle_bef << " deg" << G4endl;
    G4cout << "Angle between polDir & ptcl momDir after scattering: " << pol_mom_Angle_aft << " deg" << "\n\n";

    //G4cout << "\n\ndirectional Angle between PRF_Y and ptcl Mom: " << a
      
    //eta = acos((scatPlane.dot(prfX_before))/(scatPlane.mag()*prfX_before.mag()));
    //eta = acos((Plane.dot(prfX_before))/(scatPlane.mag()*prfX_before.mag()));
    //eta = scatPlane.angle(prfX_before);
    //eta = polPlane.angle(scatPlane);
    
    //eta = acos((scatPlane.dot(prfY_before))/(scatPlane.mag()*prfY_before.mag()));
    //eta = acos((scatPlane.dot(vin.cross(prfX_before)))/(scatPlane.mag()*(vin.cross(prfX_before)).mag()));
    //eta = scatPlane.azimAngle(prfY_before, vin);
    //eta = vout.azimAngle(prfX_before, vin);

    //polDir =

    //================ working eta===============
    //eta = std::atan2(vin.dot(prfX_before.cross(vout)), prfX_before.dot(vout));
    //===========================================

    eta = std::atan2(vin.dot(prePol.cross(vout)), prePol.dot(vout));
    //eta = std::atan2(vin.dot(polPlane.cross(vout)), polPlane.dot(vout));

    //test = G4ThreeVector()
    //eta = std::atan2(vin.dot(prfX_before.cross(vout)), prfX_before.dot(vout));
    //eta = std::atan2(vin.dot(prfY_before.cross(vout)), prfY_before.dot(vout));
    //eta = std::atan2(vin.dot(prfY_before.cross(vout)), prfY_before.dot(vout));
    
    G4ThreeVector e1 = vin.orthogonal().unit(); // E-field direction
    G4ThreeVector n_s = vin.cross(vout).unit();
    G4ThreeVector n_p = vin.cross(e1).unit();

    G4double cosPhi = n_p.dot(n_s);
    G4double sinPhi = vin.dot(e1.cross(vout)); // Gives sign

    //eta = std::atan2(sinPhi, cosPhi);  

    //eta = std::atan2(vin.dot(prePol.cross(vout)), prePol.dot(vout));

    //eta = std::atan2(vin.dot(vrandom.cross(vout)), vrandom.dot(vout));
    //eta = std::atan2(vin.dot(prfY_before.cross(vout)), prfY_before.dot(vout));

    //eta = (CLHEP::pi/2)-scatPlane.azimAngle(prfY_before, vin);

    G4cout << "vin mag: " << vin.mag() << G4endl;
    G4cout << "vout mag: " << vout.mag() << G4endl;
    G4cout << "prfX mag: " << prfX_before.mag() << G4endl;
    

    G4cout << "Angle between polPlane & ScatPlane: " << eta*180/CLHEP::pi  << " deg" << G4endl;

    //G4cout << "\n\neta is: " << eta << G4endl;
     G4cout << "\n*****" << "\n\n";

     // G4double xdir = vout * G4PolarizationHelper::GetParticleFrameX(vin);
     // G4double ydir = vout * G4PolarizationHelper::GetParticleFrameY(vin);
     
     // eta1 = std::atan2(ydir, xdir);

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
    //hit->SetEta(eta1);

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
