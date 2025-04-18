#include "run.hh"

MyRunAction::MyRunAction()
{
  // G4AnalysisManager *man = G4AnalysisManager::Instance();
  
  // man->CreateNtuple("Hits", "Hits");
  // man->CreateNtupleIColumn("fEvent");
  // man->CreateNtupleDColumn("fX");
  // man->CreateNtupleDColumn("fY");
  // man->CreateNtupleDColumn("fZ");
  // man->FinishNtuple(0);

  // man->CreateNtuple("sumX", "sumX");
  // man->CreateNtupleDColumn("sumX");
  // man->FinishNtuple(1);

  // man->CreateNtuple("position", "position");
  // man->CreateNtupleDColumn("XCoord");
  // man->FinishNtuple(2);

  // man->CreateNtuple("Scoring", "Scoring");
  // man->CreateNtupleDColumn("fEdep");
  // man->FinishNtuple(3);

  // tree = new TTree("tree","HitInfo");
  // tree->Branch("positionX", &X);
  // tree->Branch("positionY", &Y);
  // tree->Branch("positionZ", &Z);
  // tree->Branch("Energy", &E);
  // tree->Branch("Total_Edep", &Total_E);
  // tree->Branch("Edep_Compt", &Total_Compt_Edep);
  // tree->Branch("Edep_Photo", &Photo_Edep);
  // tree->Branch("diff_Edep_ComptPhoto", &diff_edep_ComptPhoto);
  // tree->Branch("nOpticalPhotons", &nOpticalPhotons);
  // tree->Branch("OptPhoOnDet", &nOptPhoOnDetEnd);
  // tree->Branch("OptPho_PosX", &Opt_Photon_PosX);
  // tree->Branch("OptPho_PosY", &Opt_Photon_PosY);
  // tree->Branch("OptPho_PosZ", &Opt_Photon_PosZ);
  // tree->Branch("OptPho_Energy", &Opt_Photon_Energy);
  // tree->Branch("OptPho_Time", &Opt_Photon_Time);
  // tree->Branch("nRefraction", &nRefraction);
  // tree->Branch("nReflection", &nReflection);
  // tree->Branch("nTIR", &nTIR);

}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  tree = new TTree("tree","HitInfo");
  tree->Branch("positionX", &X);
  tree->Branch("positionY", &Y);
  tree->Branch("positionZ", &Z);
  tree->Branch("Energy", &E);
  tree->Branch("Total_Edep", &Total_E);
  tree->Branch("Edep_Compt", &Total_Compt_Edep);
  tree->Branch("Edep_Photo", &Photo_Edep);
  tree->Branch("diff_Edep_ComptPhoto", &diff_edep_ComptPhoto);
  tree->Branch("nOpticalPhotons", &nOpticalPhotons);
  tree->Branch("OptPhoOnDet", &nOptPhoOnDetEnd);
  tree->Branch("OptPho_PosX", &Opt_Photon_PosX);
  tree->Branch("OptPho_PosY", &Opt_Photon_PosY);
  tree->Branch("OptPho_PosZ", &Opt_Photon_PosZ);
  tree->Branch("OptPho_Energy", &Opt_Photon_Energy);
  tree->Branch("OptPho_Time", &Opt_Photon_Time);
  tree->Branch("nRefraction", &nRefraction);
  tree->Branch("nReflection", &nReflection);
  tree->Branch("nTIR", &nTIR);


  G4int runID = run->GetRunID();

  std::stringstream strRunID;
  strRunID << runID;

  //man->OpenFile("output"+strRunID.str()+".root");

  hfile = hfile = TFile::Open("test.root","RECREATE");

  
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  // G4AnalysisManager *man = G4AnalysisManager::Instance();
  
  // man->Write();
  // man->CloseFile();

  tree->Write();
  //tree->Print();
    
  hfile->Close();

}

