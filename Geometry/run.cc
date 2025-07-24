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
  tree->Branch("Hit_Edep", &HitEdep);
  tree->Branch("Hit_PositionX", &HitPosX);
  tree->Branch("Hit_PositionY", &HitPosY);
  tree->Branch("Hit_PositionZ", &HitPosZ);
  tree->Branch("Hit_Time", &HitTime);
  tree->Branch("Hit_TrkLen", &HitTrklen);
  tree->Branch("Hit_DetId", &HitDetId);
  tree->Branch("Hit_ProcId", &HitProcId);
  tree->Branch("Hit_ScatAngle", &HitScatAngle);
  tree->Branch("Hit_Eta", &HitEta);
  tree->Branch("Hit_Pol0", &HitPol0);
  tree->Branch("Hit_Pol1", &HitPol1);
  tree->Branch("Hit_Pol2", &HitPol2);
  tree->Branch("Hit_Eout", &HitEout);
  tree->Branch("Hit_Ein", &HitEin);

  hfile = hfile = TFile::Open("test.root","RECREATE");
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();
  //tree->Print();
    
  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;
}

