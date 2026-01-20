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

  tree->Branch("nOpticalPhotons", &nOpPhotons);
  tree->Branch("Hit_Edep", &HitEdep);
  tree->Branch("Hit_PositionX", &HitPosX);
  tree->Branch("Hit_PositionY", &HitPosY);
  tree->Branch("Hit_PositionZ", &HitPosZ);
  tree->Branch("Hit_Time", &HitTime);
  tree->Branch("Hit_DetId", &HitDetId);
  tree->Branch("Hit_GunId", &HitGunId);
  tree->Branch("Hit_ProcId", &HitProcId);
  tree->Branch("Hit_ScatAngle", &HitScatAngle);
  tree->Branch("Hit_Eta", &HitEta);
  tree->Branch("Hit_Eout", &HitEout);
  tree->Branch("Hit_Ein", &HitEin);
  tree->Branch("right_module", Module1, "Module1[8][8][7]/F");
  tree->Branch("left_module", Module2, "Module2[8][8][7]/F");

  tree->Branch("Edep_truch", &Edep_truth);
  tree->Branch("nPtcl_OutDet", &nPtcl_out);
  tree->Branch("E_outPtcl", &E_outPtcl);
  tree->Branch("E_verify", &E_verify);
  tree->Branch("nSec_e", &nSec_e);
  tree->Branch("nSec_pho", &nSec_pho);
  tree->Branch("E_elec_out", &E_eOut);
  tree->Branch("E_pho_out", &E_phoOut);
  tree->Branch("ScatTheta", &scat_theta);

  tree->Branch("rel_PolAngle", &relPol_Angle);
  
  hfile = hfile = TFile::Open("Uncorrelated_test_4.5M_2Mod_2Gun.root","RECREATE");
}



void MyRunAction::EndOfRunAction(const G4Run*)
{
  tree->Write();
  //tree->Print();
    
  hfile->Close();
  G4cout << "Run Ended !!" << G4endl;

  // G4cout << "\v============= MATRIX 1 (GenRun): =============" << G4endl;
  // for (int i = 0; i < 8; i++) {          // Loop over rows
  //   for (int j = 0; j < 8; j++) {      // Loop over columns
  //     G4cout << std::setw(6) << Module1[i][j][nGenOp] << " ";
  //   }
  //   G4cout << G4endl;                      // New line after each row
  // }
  
 //  G4cout << "\n\n\v============= MATRIX 2: =============" << G4endl;
//   for (int i = 0; i < 8; i++) {          // Loop over rows
//     for (int j = 0; j < 8; j++) {      // Loop over columns
//       G4cout << std::setw(6) << Module2[i][j] << " ";
//     }
//     G4cout << G4endl;                      // New line after each row
//   }  
}
