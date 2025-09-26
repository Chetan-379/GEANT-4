#include "run.hh"

MyRunAction::MyRunAction()
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();  
  man->CreateNtuple("OpPho_Ana", "OpPho_Ana");
  man->CreateNtupleDColumn("OpPho_Lambda");
  man->CreateNtupleDColumn("OpPho_Z");
  man->FinishNtuple(0);

  man->CreateNtuple("OpPho_End", "OpPho_End");
  man->CreateNtupleDColumn("OpPho_Lambda");
  man->CreateNtupleDColumn("CreatPos_Z");
  man->CreateNtupleDColumn("OpPho_PosZ");
  man->FinishNtuple(1);

  man->CreateNtuple("EvtInfo", "EvtInfo");
  man->CreateNtupleIColumn("nOpPho");
  man->CreateNtupleIColumn("nOpPho_End");
  man->FinishNtuple(2);

  


  // tree = new TTree("tree","OpticalPho_Ana");
  // tree->Branch("OpPho_Lambda", &OpPho_lmbda);
  // tree->Branch("OpPho_CreateZ", &OpPho_CreateZ);

  // hfile = hfile = TFile::Open("test.root","RECREATE");
}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "run Started!!" << G4endl;
  
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  man->OpenFile("test.root");
 
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();  
  man->Write();
  man->CloseFile();
  


  // tree->Write();
  // hfile->Close();
  
  G4cout << "Run Ended !!" << G4endl;
}
