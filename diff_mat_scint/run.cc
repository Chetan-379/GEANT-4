#include "run.hh"

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "run Started!!" << G4endl;
  
  tree = new TTree("tree","HitInfo");
  tree->Branch("positionX", &X);
  tree->Branch("positionY", &Y);
  tree->Branch("positionZ", &Z);
  tree->Branch("Energy", &E);
  tree->Branch("Total_Edep", &Total_E);
  tree->Branch("Edep_Compt", &Total_Compt_Edep);
  tree->Branch("Edep_Photo", &Photo_Edep);
  tree->Branch("nOpticalPhotons", &nOpticalPhotons);
  tree->Branch("OptPhoOnDet", &nOptPhoOnDetEnd);
  tree->Branch("OptPho_PosX", &Opt_Photon_PosX);
  tree->Branch("OptPho_PosY", &Opt_Photon_PosY);
  tree->Branch("OptPho_PosZ", &Opt_Photon_PosZ);
  tree->Branch("OptPho_Energy", &Opt_Photon_Energy);
  tree->Branch("OptPho_Time", &Opt_Photon_Time);
  tree->Branch("nCompton", &nCompton);
  tree->Branch("theta", &theta);

  // const PrimaryGeneratorAction* genAction =
  //   static_cast<const PrimaryGeneratorAction*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  //  G4double energy = genAction->GetGunEnergy();
  //  G4cout << energy << G4endl;

  // G4int RunID = run->GetRunID();
  // G4cout << RunID << G4endl;
  // G4String depth, fileName;

  // if (RunID ==0) depth = "1cm";
  // if (RunID ==1) depth = "2.5cm";
  // if (RunID ==2) depth = "5cm";
  // if (RunID ==3) depth = "10cm";

  //fileName = "test_" + depth + G4String(energy)+"keV.root";
  // fileName = "test_" + depth + ".root";
  
  hfile = TFile::Open("test.root","RECREATE");
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();
  //tree->Print();
    
  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;
}

