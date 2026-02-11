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
  tree->Branch("nOpGAGG_quad", nOpG_quad, "nOpG_quad[4]/I");
  tree->Branch("nOpLYSO_quad", nOpL_quad, "nOpL_quad[4]/I");
  
  tree->Branch("RandomSeed", &gRandomSeed, "RandomSeed/L");

  //h_Op_lmbda = new TH1F("Op_lambda","Op_lambda",1000,0,1000);

  h_zVslmbda = new TH2F("Op_lambda","Op_lambda", 100,15,35, 1000,0,1000);
  h_zVslmbda->GetXaxis()->SetTitle("CreateZ (mm)");
  h_zVslmbda->GetYaxis()->SetTitle("Op_lambda (nm)");
  //output root file------------------
  if(gOutputFileName != " ") hfile = hfile = TFile::Open((gOutputFileName).c_str(),"RECREATE");
  else hfile = hfile = TFile::Open("test_vis.root","RECREATE");

  G4cout << "----------------------Seed: " << gRandomSeed << "---------------" << G4endl;
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();
  //h_Op_lmbda->Write();
  h_zVslmbda->Write();
  

  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;
}
