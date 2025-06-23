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
  tree->Branch("Hit_Edep", &InHitEdep);
  tree->Branch("Hit_PositionX", &InHitPosX);
  tree->Branch("Hit_PositionY", &InHitPosY);
  tree->Branch("Hit_PositionZ", &InHitPosZ);
  tree->Branch("Hit_Time", &InHitTime);
  tree->Branch("Hit_TrkLen", &InHitTrklen);
  tree->Branch("Hit_DetId", &InHitDetId);

  hfile = hfile = TFile::Open("test.root","RECREATE");
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();
  //tree->Print();
    
  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;
}

