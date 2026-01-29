#include "run.hh"

extern long gRandomSeed;
extern std::string gOutputFileName;

MyRunAction::MyRunAction()
{}

MyRunAction::~MyRunAction()
{}

void MyRunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "run Started!!" << G4endl;
  
  tree = new TTree("tree","HitInfo");
  tree->Branch("nOp_GAGG", &nOpGAGG);
  tree->Branch("nOp_LYSO", &nOpLYSO);
  tree->Branch("nOp_GAGG_truth", &nOpGAGG_truth);
  tree->Branch("nOp_LYSO_truth", &nOpLYSO_truth);
  tree->Branch("Edep_truth", &Edep_truth);
  tree->Branch("Edep_GAGG_truth", &Edep_G_truth);
  tree->Branch("Edep_LYSO_truth", &Edep_L_truth);
  
  tree->Branch("RandomSeed", &gRandomSeed, "RandomSeed/L");
  
  //output root file------------------
  if(gOutputFileName != " ") hfile = hfile = TFile::Open((gOutputFileName).c_str(),"RECREATE");
  else hfile = hfile = TFile::Open("test_vis.root","RECREATE");

  G4cout << "----------------------Seed: " << gRandomSeed << "---------------" << G4endl;
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();

  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;
}
