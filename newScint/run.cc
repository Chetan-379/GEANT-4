#include "run.hh"

MyRunAction::MyRunAction()
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  
  man->CreateNtuple("OpPho_Ana", "OpPho_Ana");
  man->CreateNtupleDColumn("OpPho_Lambda");
  man->CreateNtupleDColumn("OpPho_Z");

  man->FinishNtuple(0);

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

  G4cout << "Run Ended !!" << G4endl;
}
