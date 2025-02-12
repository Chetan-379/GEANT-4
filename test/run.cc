#include "run.hh"

MyRunAction::MyRunAction()
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  
  man->CreateNtuple("Hits", "Hits");
  man->CreateNtupleIColumn("fEvent");
  man->CreateNtupleDColumn("fX");
  man->CreateNtupleDColumn("fY");
  man->CreateNtupleDColumn("fZ");
  man->FinishNtuple(0);

  man->CreateNtuple("sumX", "sumX");
  man->CreateNtupleDColumn("sumX");
  man->FinishNtuple(1);

  man->CreateNtuple("position", "position");
  man->CreateNtupleDColumn("XCoord");
  man->FinishNtuple(2);

  //std::cout << man << std::endl;

  tree = new TTree("tree","HitInfo");
  tree->Branch("positionX", &X);
  tree->Branch("positionY", &Y);
  tree->Branch("positionZ", &Z);
  tree->Branch("Energy", &E);
  tree->Branch("Total_Edep", &Total_E);
  //x=p;
  //tree = nullptr;
  //Xcoord = 0;
  

}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  G4int runID = run->GetRunID();

  std::stringstream strRunID;
  strRunID << runID;

  man->OpenFile("output"+strRunID.str()+".root");

  //TFile *hfile;
  //hfile = TFile::Open("position.root","RECREATE");
  //TFile *hfile;
  hfile = hfile = TFile::Open("test.root","RECREATE");
  
  
  //tree->Branch("positionX", &Xcoord, "positionX/I");
  
 
  //std::cout << man << std::endl;
}

//void MyRunAction::FillTree(std::vector<double> Xvec){
 
//}


// void MyRunAction::FillTreee(){
//   tree->Fill(); 
//}


void MyRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  
  man->Write();
  man->CloseFile();

  //TTree *tree;
  //tree->Write();
  //tree->Close();
  //  tree->Branch("positionX", &x, "positionX/I");
  //tree->Fill();
  tree->Write();
  tree->Print();
    
  hfile->Close();
  //std::cout << man << std::endl;
}

