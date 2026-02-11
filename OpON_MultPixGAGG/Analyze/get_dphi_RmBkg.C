void get_dphi_RmBkg() {
    // TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_5MEvts_correct_DetPos_out.root");
    // TFile *f2 = TFile::Open("out_root_files/OpOFF_Uncorrelated_5MEvts_correct_DetPos_out.root");

    //TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_theta_70_90_5MEvts_correct_DetPos_out.root");
  // TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_theta_70_90_5MEvts_correct_DetPos_correct_pol_out.root");
  // TFile *f2 = TFile::Open("out_root_files/OpOFF_Uncorrelated_theta_70_90_5MEvts_correct_DetPos_out.root");

  // TFile *f1 = TFile::Open("out_root_files/plotting_reco_info_Correlated_OpON_4.5MEvts.root");
  // TFile *f2 = TFile::Open("out_root_files/Uncorrelated_OpON_4.5MEvts.root");

 TFile *f1 = TFile::Open("check_correlated_50MEvts.root");
 TFile *f2 = TFile::Open("check_Uncorrelated_50MEvts.root");

  // TFile *f1 = TFile::Open("check_correlated.root");
  // TFile *f2 = TFile::Open("check_Uncorrelated.root");


 
  // TH1F *h1 = (TH1F*) f1->Get("dPhi_truth_70_90");
  // TH1F *h2 = (TH1F*) f2->Get("dPhi_truth_70_90");

 // TH1F *h1 = (TH1F*) f1->Get("dPhi_reco_70_90");
 // TH1F *h2 = (TH1F*) f2->Get("dPhi_reco_70_90");
 
  TH1F *h1 = (TH1F*) f1->Get("dPhi_truth_inc_dexc");
  TH1F *h2 = (TH1F*) f2->Get("dPhi_truth_inc_dexc");

  // TH1F *h1 = (TH1F*) f1->Get("dPhi_truth_70_90_dexc");
  // TH1F *h2 = (TH1F*) f2->Get("dPhi_truth_70_90_dexc");

  // h1->Scale(1.0 / (h1->Integral()));  
  // h2->Scale(1.0 / (h2->Integral()));
 
  h1->Sumw2();
  h2->Sumw2();
  
  h1->Rebin(5);
  h2->Rebin(5);
  
  TH1F *hDiv = (TH1F*) h1->Clone("hDiv");
  hDiv->Divide(h2);
  
  //hDiv->Rebin(5);
  
  TFile *fOut = TFile::Open("divided.root", "RECREATE");
  h1->Write();
  h2->Write();
  hDiv->Write();
  
  fOut->Close();
  
  f1->Close();
  f2->Close();      
}
