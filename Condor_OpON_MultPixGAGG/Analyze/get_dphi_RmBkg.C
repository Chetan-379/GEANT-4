void get_dphi_RmBkg() {
    // TFile *f1 = TFile::Open("out_root_files/correlated_theta_70_90_5MEvts.root");
    // TFile *f2 = TFile::Open("out_root_files/Uncorrelated_theta_70_90_5MEvts.root");

    // TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_5MEvts_correct_DetPos_out.root");
    // TFile *f2 = TFile::Open("out_root_files/OpOFF_Uncorrelated_5MEvts_correct_DetPos_out.root");

    //TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_theta_70_90_5MEvts_correct_DetPos_out.root");
  TFile *f1 = TFile::Open("out_root_files/OpOFF_Correlated_theta_70_90_5MEvts_correct_DetPos_correct_pol_out.root");
  TFile *f2 = TFile::Open("out_root_files/OpOFF_Uncorrelated_theta_70_90_5MEvts_correct_DetPos_out.root");
  
  TH1F *h1 = (TH1F*) f1->Get("diff_Phi_truth");
  TH1F *h2 = (TH1F*) f2->Get("diff_Phi_truth");
  
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
