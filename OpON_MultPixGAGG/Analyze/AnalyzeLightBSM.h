#ifndef AnalyzeLightBSM_H
#define AnalyzeLightBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include <TProfile.h>
#include "TFile.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include "TGaxis.h"


#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;

class AnalyzeLightBSM : public NtupleVariables{

public:
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *detType="PbWO4");
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *);
  void     BookHistogram(const char *);
  void FillHistogram(double, double, double, double, double, double, double, double);
  string MatName(int);
  void print(Long64_t);
  float GetScatAngle(float);
  float GetScatElecE(float);
  float Deg(float);
  float Rad(float);

  float Pi = TMath::Pi();
  //define histograms here
  TH1D *h_Compt_Edep, *h_Photo_Edep, *h_Total_Edep, *h_Total_Edep_fine_binned;
  TH2D *h_ComptVsPhoto_Edep;
  TH2D *h_ScatAng_SDvsAna;
  TH1D *h_diff_Ana_SD;
  TH1I *h_nHits;
  
  TH1D *h_theta, *h_eta;
  TH1D *h_hitTime;
  TH1F *h_hitPosX, *h_hitPosY, *h_hitPosZ;
  TH1F *h_Eout, *h_Ein;
  TH2D *h_theta_Eout;

  TH1I *h_nOpPho_crys_inc[2], *h_nOpPho_gen_crys_inc[2], *h_nOpPho_mdl_inc[2], *h_nOpPho_gen_mdl_inc[2];
  TH1I *h_nOpPho_crys_indv[2], *h_nOpPho_gen_crys_indv[2], *h_nOpPho_mdl_indv[2], *h_nOpPho_gen_mdl_indv[2];
  TH1I *h_nOpPho_crys_sum[2], *h_nOpPho_gen_crys_sum[2], *h_nOpPho_mdl_sum[2], *h_nOpPho_gen_mdl_sum[2];
  

  TH1F *h_DetX_crys[2], *h_Edep_truth_crys[2], *h_Edep_truth_mdl[2], *h_Edep_reco_mdl[2], *h_Edep_reco_crys[2];

  TH2D *h_nOp_vs_Edep_crys[2], *h_nOp_vs_Edep_mdl[2];

  TH1F* h_Edep_diff[2], *h_Edep_scat_truth[2], *h_Edep_abs_truth[2], *h_Edep_scat_reco[2], *h_Edep_abs_reco[2];

  TH2F* h_Edep_true_Vs_reco_crys[2], *h_Edep_scatVsabs_truth[2], *h_Edep_scatVsabs_reco[2], *h_DetXvsDetY_truth[2], *h_DetXvsDetY_reco[2];

  TH1I *h_nHits_truth_mdl[2], *h_nHits_reco_mdl[2];

  TH1F* h_theta_reco[2], *h_theta_truth[2], *h_scatElec_E[2], *h_phi_truth[2], *h_phi_truth_EvtMix[2], *h_phi_reco[2], *h_dPix_truth[2], *h_dPix_reco[2];

  TH2F* h_elecE_vs_Edep_crys[2], *h_elecE_vs_scatEdep[2];

  TH1F *h_dPhi_truth_inc, *h_dPhi_truth_inc_dexc, *h_dPhi_reco_inc, *h_dPhi_truth_EvtMix_inc, *h_dPhi_truth_corrected_inc;
  TH1F *h_dPhi_truth_70_90, *h_dPhi_truth_70_90_dexc, *h_dPhi_reco_70_90, *h_dPhi_truth_EvtMix_70_90, *h_dPhi_truth_corrected_70_90;
  
  TH2F *h_Phi1_vs_Phi2_truth, *h_Phi1_vs_Phi2_reco;

  TH2F *h_test2d;
  TH1F *h_Edep_sum_reco;
   
  TFile *oFile;

  enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep = 4, E_elec=5, tryVar =6};
};

enum ModIdx{r=0, l=1};

enum sigma_check{inc=0, hitIndv=1, hitSum=2};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //TGaxis::SetMaxDigits(2);

  h_Compt_Edep = new TH1D("Compton_Edep","Edep_via_Compton",100,0,1.1);
  h_Photo_Edep = new TH1D("Photo_Edep","Edep_via_Photoelecric",100,0,1.1);
  h_Total_Edep = new TH1D("Total_Edep","Total_Edep",100,0,1.1);
   
  h_ComptVsPhoto_Edep = new TH2D("ComptVsPho","Compt_vs_Photo_Edep",100,0.,1.1, 100,0.,1.1);
  h_ComptVsPhoto_Edep->SetXTitle("Comp_Edep");
  h_ComptVsPhoto_Edep->SetYTitle("PhotoElectric_Edep");

  double bin_edges[16];
  for (int j = 0; j <= 15; ++j) {
    bin_edges[j] = j - 0.5;
  }    
  h_nHits = new TH1I("nHits", "nHits", 15, bin_edges);

  h_hitTime = new TH1D("h_hitTime", "h_hitTime", 1000,0 ,10);
  
  h_hitPosX = new TH1F("h_hitPosX", "h_hitPosX", 250,-700 ,700);
  h_hitPosY = new TH1F("h_hitPosY", "h_hitPosY", 250,-700 ,700);
  h_hitPosZ = new TH1F("h_hitPosZ", "h_hitPosZ", 250,-700 ,700);      
  
  h_theta = new TH1D("h_theta", "h_theta", 100,0, 200);
  h_theta->GetXaxis()->SetTitle("theta (deg)");
  
  h_eta = new TH1D("h_eta", "h_eta", 200, -200, 200);
  h_eta->GetXaxis()->SetTitle("eta (deg)");

  h_Eout = new TH1F("h_Eout", "h_Eout", 260, 0, 0.520);
  h_Eout->GetXaxis()->SetTitle("Eout (MeV)");
  
  h_Ein = new TH1F("h_Ein", "h_Ein", 260, 0, 0.520);
  h_Ein->GetXaxis()->SetTitle("Ein (MeV)");

  //------------Histograms for the crystals and modules----------------------
  vector<string> ModName = {"right", "left"};
  
  char hname_nOpPho_cryst_inc[50], hname_nOpPho_gen_cryst_inc[50], hname_nOpPho_mdl_inc[50], hname_nOpPho_gen_mdl_inc[50], hname_nOpPho_cryst_indv[50], hname_nOpPho_gen_cryst_indv[50], hname_nOpPho_mdl_indv[50], hname_nOpPho_gen_mdl_indv[50], hname_nOpPho_cryst_sum[50], hname_nOpPho_gen_cryst_sum[50], hname_nOpPho_mdl_sum[50], hname_nOpPho_gen_mdl_sum[50];

  char hname_DetPosX_crys[50], hname_Edep_truth_cryst[50], hname_Edep_reco_cryst[50], hname_nOp_vs_Edep_cryst[50], hname_Edep_truth_mdl[50], hname_Edep_reco_mdl[50], hname_nOp_vs_Edep_mdl[50], hname_Edep_diff[50], hname_Edep_TrueVsReco_cryst[50], hname_nHits_truth_mdl[50], hname_nHits_reco_mdl[50], hname_DetXvsDetY_truth[50], hname_DetXvsDetY_reco[50], hname_theta_reco[50], hname_theta_truth[50], hname_scat_ElecVsEdep[50], hname_scat_ElecVSscat_Edep[50], hname_Edep_scat_truth[50], hname_Edep_abs_truth[50],  hname_Edep_scat_reco[50], hname_Edep_abs_reco[50], hname_Edep_scatVsabs_truth[50], hname_Edep_scatVsabs_reco[50], hname_scat_elec_E[50], hname_phi_truth[50], hname_phi_truth_EvtMix[50], hname_phi_reco[50], hname_dPix_truth[50], hname_dPix_reco[50];

  for (int i=0; i<2; i++){
  
  sprintf(hname_nOpPho_cryst_inc, "%s_nOpPho_cryst_inc",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_cryst_inc, "%s_nOpPho_gen_cryst_inc",ModName[i].c_str());
  sprintf(hname_nOpPho_mdl_inc, "%s_nOpPho_mdl_inc",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_mdl_inc, "%s_nOpPho_gen_mdl_inc",ModName[i].c_str());
  
  sprintf(hname_nOpPho_cryst_indv, "%s_nOpPho_cryst_indv",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_cryst_indv, "%s_nOpPho_gen_cryst_indv",ModName[i].c_str());
  sprintf(hname_nOpPho_mdl_indv, "%s_nOpPho_mdl_indv",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_mdl_indv, "%s_nOpPho_gen_mdl_indv",ModName[i].c_str());

  sprintf(hname_nOpPho_cryst_sum, "%s_nOpPho_cryst_sum",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_cryst_sum, "%s_nOpPho_gen_cryst_sum",ModName[i].c_str());
  sprintf(hname_nOpPho_mdl_sum, "%s_nOpPho_mdl_sum",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_mdl_sum, "%s_nOpPho_gen_mdl_sum",ModName[i].c_str());
  

  sprintf(hname_DetPosX_crys, "%s_DetPosX_crys", ModName[i].c_str());
  sprintf(hname_Edep_truth_cryst, "%s_Edep_truth_cryst",ModName[i].c_str());
  sprintf(hname_Edep_reco_cryst, "%s_Edep_reco_cryst",ModName[i].c_str());

  sprintf(hname_Edep_truth_mdl, "%s_Edep_truth_mdl",ModName[i].c_str());
  sprintf(hname_Edep_reco_mdl, "%s_Edep_reco_mdl",ModName[i].c_str());

  sprintf(hname_nHits_truth_mdl, "%s_nHits_truth_mdl",ModName[i].c_str());
  sprintf(hname_nHits_reco_mdl, "%s_nHits_reco_mdl",ModName[i].c_str());

  sprintf(hname_dPix_truth, "%s_pixel_distance_truth",ModName[i].c_str());
  sprintf(hname_dPix_reco, "%s_pixel_distance_reco",ModName[i].c_str());

  sprintf(hname_Edep_scat_truth, "%s_Edep_scat_truth",ModName[i].c_str());
  sprintf(hname_Edep_abs_truth, "%s_Edep_abs_truth",ModName[i].c_str());
  sprintf(hname_Edep_scatVsabs_truth, "%s_Edep_scatVsabs_truth",ModName[i].c_str());

  sprintf(hname_Edep_scat_reco, "%s_Edep_scat_reco",ModName[i].c_str());
  sprintf(hname_Edep_abs_reco, "%s_Edep_abs_reco",ModName[i].c_str());
  sprintf(hname_Edep_scatVsabs_reco, "%s_Edep_reco",ModName[i].c_str());

  sprintf(hname_DetXvsDetY_truth, "%s_DetXvsDetY_truth",ModName[i].c_str());
  sprintf(hname_DetXvsDetY_reco, "%s_DetXvsDetY_reco",ModName[i].c_str());


  sprintf(hname_phi_truth, "%s_phi_truth",ModName[i].c_str());
  sprintf(hname_phi_truth_EvtMix, "%s_phi_truth_EvtMix",ModName[i].c_str());
  
  sprintf(hname_phi_reco, "%s_phi_reco",ModName[i].c_str());

  
  sprintf(hname_nOp_vs_Edep_cryst, "%s_nOp_vs_Edep_cryst",ModName[i].c_str());
  sprintf(hname_nOp_vs_Edep_mdl, "%s_nOp_vs_Edep_mdl",ModName[i].c_str());
  sprintf(hname_Edep_diff, "%s_Edep_diff",ModName[i].c_str());
  sprintf(hname_Edep_TrueVsReco_cryst, "%s_Edep_TrueVsReco_cryst",ModName[i].c_str());    
  sprintf(hname_theta_reco, "%s_theta_reco",ModName[i].c_str());
  sprintf(hname_theta_truth, "%s_theta_truth",ModName[i].c_str());
  sprintf(hname_scat_ElecVsEdep, "%s_scat_ElecVsEdep",ModName[i].c_str());
  sprintf(hname_scat_ElecVSscat_Edep, "%s_scat_ElecVSscat_Edep",ModName[i].c_str());
  
  sprintf(hname_scat_elec_E, "%s_scat_elec_E",ModName[i].c_str());
  
  
  

  
  //inc-------
  h_nOpPho_crys_inc[i] = new TH1I(hname_nOpPho_cryst_inc, hname_nOpPho_cryst_inc,3000,0,30000);
  h_nOpPho_crys_inc[i]->GetXaxis()->SetTitle("nOp_end");

  h_nOpPho_gen_crys_inc[i] = new TH1I(hname_nOpPho_gen_cryst_inc, hname_nOpPho_gen_cryst_inc,3000,0,30000);
  h_nOpPho_gen_crys_inc[i]->GetXaxis()->SetTitle("nOp_gen");
  
  h_nOpPho_mdl_inc[i] = new TH1I(hname_nOpPho_mdl_inc, hname_nOpPho_mdl_inc,3000,0,30000);
  h_nOpPho_mdl_inc[i]->GetXaxis()->SetTitle("nOp_end");
  
  h_nOpPho_gen_mdl_inc[i] = new TH1I(hname_nOpPho_gen_mdl_inc, hname_nOpPho_gen_mdl_inc,3000,0,30000);
  h_nOpPho_gen_crys_inc[i]->GetXaxis()->SetTitle("nOp_gen");

  //indv-------
  h_nOpPho_crys_indv[i] = new TH1I(hname_nOpPho_cryst_indv, hname_nOpPho_cryst_indv,3000,0,30000);
  h_nOpPho_crys_indv[i]->GetXaxis()->SetTitle("nOp_end");

  h_nOpPho_gen_crys_indv[i] = new TH1I(hname_nOpPho_gen_cryst_indv, hname_nOpPho_gen_cryst_indv,3000,0,30000);
  h_nOpPho_gen_crys_indv[i]->GetXaxis()->SetTitle("nOp_gen");
  
  h_nOpPho_mdl_indv[i] = new TH1I(hname_nOpPho_mdl_indv, hname_nOpPho_mdl_indv,3000,0,30000);
  h_nOpPho_mdl_indv[i]->GetXaxis()->SetTitle("nOp_end");
  
  h_nOpPho_gen_mdl_indv[i] = new TH1I(hname_nOpPho_gen_mdl_indv, hname_nOpPho_gen_mdl_indv,3000,0,30000);
  h_nOpPho_gen_crys_indv[i]->GetXaxis()->SetTitle("nOp_gen");

  //sum-------
  h_nOpPho_crys_sum[i] = new TH1I(hname_nOpPho_cryst_sum, hname_nOpPho_cryst_sum,3000,0,30000);
  h_nOpPho_crys_sum[i]->GetXaxis()->SetTitle("nOp_end");

  h_nOpPho_gen_crys_sum[i] = new TH1I(hname_nOpPho_gen_cryst_sum, hname_nOpPho_gen_cryst_sum,3000,0,30000);
  h_nOpPho_gen_crys_sum[i]->GetXaxis()->SetTitle("nOp_gen");
  
  h_nOpPho_mdl_sum[i] = new TH1I(hname_nOpPho_mdl_sum, hname_nOpPho_mdl_sum,3000,0,30000);
  h_nOpPho_mdl_sum[i]->GetXaxis()->SetTitle("nOp_end");
  
  h_nOpPho_gen_mdl_sum[i] = new TH1I(hname_nOpPho_gen_mdl_sum, hname_nOpPho_gen_mdl_sum,3000,0,30000);
  h_nOpPho_gen_crys_sum[i]->GetXaxis()->SetTitle("nOp_gen");

  //-----------------------------------------------------
  h_DetX_crys[i] = new TH1F(hname_DetPosX_crys, hname_DetPosX_crys,1000,-40,40);  
  h_Edep_truth_crys[i] = new TH1F (hname_Edep_truth_cryst, hname_Edep_truth_cryst, 100,0,1.1);
  h_Edep_truth_crys[i]->GetXaxis()->SetTitle("Truth Edep (MeV)");
  
  h_Edep_reco_crys[i] = new TH1F (hname_Edep_reco_cryst, hname_Edep_reco_cryst, 500,0,1.1);
  h_Edep_reco_crys[i]->GetXaxis()->SetTitle("Reco Edep (MeV)");
  
  h_nOp_vs_Edep_crys[i] = new TH2D(hname_nOp_vs_Edep_cryst, hname_nOp_vs_Edep_cryst, 100,0.,1.1, 3000,0,30000);
  h_nOp_vs_Edep_crys[i]->GetXaxis()->SetTitle("True Edep (MeV)");
  h_nOp_vs_Edep_crys[i]->GetYaxis()->SetTitle("nOpPho_end");

  h_Edep_truth_mdl[i] = new TH1F (hname_Edep_truth_mdl, hname_Edep_truth_mdl, 100,0,1.1);
  h_Edep_truth_mdl[i]->GetXaxis()->SetTitle("True Edep (MeV)");
  
  h_Edep_reco_mdl[i] = new TH1F (hname_Edep_reco_mdl, hname_Edep_reco_mdl, 100,0,1.1);
  h_Edep_reco_mdl[i]->GetXaxis()->SetTitle("Reco Edep (MeV)");
  
  h_nOp_vs_Edep_mdl[i] = new TH2D(hname_nOp_vs_Edep_mdl, hname_nOp_vs_Edep_mdl, 100,0.,1.1, 3000,0,30000);
  h_nOp_vs_Edep_mdl[i]->GetXaxis()->SetTitle("True Edep (MeV)");
  h_nOp_vs_Edep_mdl[i]->GetYaxis()->SetTitle("nOpPho_end");

  h_Edep_diff[i] = new TH1F (hname_Edep_diff, hname_Edep_diff, 1000,0.,0.5);

  h_Edep_true_Vs_reco_crys[i] = new TH2F(hname_Edep_TrueVsReco_cryst, hname_Edep_TrueVsReco_cryst, 100,0.,1.1, 100,0,1.1);
  h_Edep_true_Vs_reco_crys[i]->GetXaxis()->SetTitle("Edep_{truth} (MeV)");
  h_Edep_true_Vs_reco_crys[i]->GetYaxis()->SetTitle("Edep_{reco} (MeV)");

  h_nHits_truth_mdl[i] = new TH1I(hname_nHits_truth_mdl, hname_nHits_truth_mdl, 10,0,10);
  h_nHits_reco_mdl[i] = new TH1I(hname_nHits_reco_mdl, hname_nHits_reco_mdl, 10,0,10);

  h_DetXvsDetY_truth[i] = new TH2F(hname_DetXvsDetY_truth, hname_DetXvsDetY_truth, 100,-40,40, 100,-40,40);
  h_DetXvsDetY_truth[i]->GetXaxis()->SetTitle("DetX_truth (mm)");
  h_DetXvsDetY_truth[i]->GetYaxis()->SetTitle("DetY_truth (mm)");

  h_DetXvsDetY_reco[i] = new TH2F(hname_DetXvsDetY_reco, hname_DetXvsDetY_reco, 100,-40,40, 100,-40,40);
  h_DetXvsDetY_reco[i]->GetXaxis()->SetTitle("DetX_reco (mm)");
  h_DetXvsDetY_reco[i]->GetYaxis()->SetTitle("DetY_reco (mm)");


  h_theta_reco[i] = new TH1F(hname_theta_reco, hname_theta_reco, 400,0,200);
  h_theta_truth[i] = new TH1F(hname_theta_truth, hname_theta_truth, 400,0,200);

  h_elecE_vs_Edep_crys[i] = new TH2F(hname_scat_ElecVsEdep, hname_scat_ElecVsEdep, 100,0.,1.1, 100,0.,1.1);
  h_elecE_vs_Edep_crys[i]->GetXaxis()->SetTitle("Edep_truth (MeV)");
  h_elecE_vs_Edep_crys[i]->GetYaxis()->SetTitle("E_scat_elec (MeV)");

  h_elecE_vs_scatEdep[i] = new TH2F(hname_scat_ElecVSscat_Edep, hname_scat_ElecVSscat_Edep, 100,0.,1.1, 100,0.,1.1);
  h_elecE_vs_scatEdep[i]->GetXaxis()->SetTitle("Edep_scat (MeV)");
  h_elecE_vs_scatEdep[i]->GetYaxis()->SetTitle("E_scat_elec (MeV)");


  h_Edep_scat_truth[i] = new TH1F(hname_Edep_scat_truth, hname_Edep_scat_truth, 100,0,0.8);
  h_Edep_abs_truth[i] = new TH1F(hname_Edep_abs_truth, hname_Edep_abs_truth, 100,0,0.8);

  h_Edep_scat_reco[i] = new TH1F(hname_Edep_scat_reco, hname_Edep_scat_reco, 100,0,0.8);
  h_Edep_abs_reco[i] = new TH1F(hname_Edep_abs_reco, hname_Edep_abs_reco, 100,0,0.8);


  h_Edep_scatVsabs_truth[i] = new TH2F(hname_Edep_scatVsabs_truth, hname_Edep_scatVsabs_truth, 100,0,0.8, 100,0,0.8);
  h_Edep_scatVsabs_truth[i]->GetXaxis()->SetTitle("scat_Edep_truth (MeV)");
  h_Edep_scatVsabs_truth[i]->GetYaxis()->SetTitle("abs_Edep_truth (MeV)");

    h_Edep_scatVsabs_reco[i] = new TH2F(hname_Edep_scatVsabs_reco, hname_Edep_scatVsabs_reco, 100,0,0.8, 100,0,0.8);
  h_Edep_scatVsabs_reco[i]->GetXaxis()->SetTitle("scat_Edep_reco (MeV)");
  h_Edep_scatVsabs_reco[i]->GetYaxis()->SetTitle("abs_Edep_reco (MeV)");

  h_scatElec_E[i]=new TH1F(hname_scat_elec_E, hname_scat_elec_E, 100,0,0.8);

  h_phi_truth[i] = new TH1F(hname_phi_truth, hname_phi_truth, 100,-200,200);
  h_phi_truth_EvtMix[i] = new TH1F(hname_phi_truth_EvtMix, hname_phi_truth_EvtMix, 100,-200,200);

  h_phi_reco[i] = new TH1F(hname_phi_reco, hname_phi_reco, 100,-200,200);

  h_dPix_truth[i] = new TH1F(hname_dPix_truth, hname_dPix_truth, 40,0,40);
  h_dPix_reco[i] = new TH1F(hname_dPix_reco, hname_dPix_reco, 40,0,40);
  }
  //---------------------------
  
  
  
  h_Phi1_vs_Phi2_truth = new TH2F("Phi1_vs_Phi2_truth", "Phi1_vs_Phi2_truth", 100,-200,200, 100,-200,200);
  h_Phi1_vs_Phi2_truth->GetXaxis()->SetTitle("Phi1_truth (deg)");
  h_Phi1_vs_Phi2_truth->GetYaxis()->SetTitle("Phi2_truth (deg)");

  h_Phi1_vs_Phi2_reco = new TH2F("Phi1_vs_Phi2_reco", "Phi1_vs_Phi2_reco", 100,-200,200, 100,-200,200);
  h_Phi1_vs_Phi2_reco->GetXaxis()->SetTitle("Phi1_reco (deg)");
  h_Phi1_vs_Phi2_reco->GetYaxis()->SetTitle("Phi2_reco (deg)");

  //------------distributions for dPhi
  h_dPhi_truth_inc = new TH1F("dPhi_truth_inc", "dPhi_truth_inc", 100,-200,200);
  h_dPhi_truth_inc->GetXaxis()->SetTitle("dPhi_truth (deg)");

  h_dPhi_truth_inc_dexc = new TH1F("dPhi_truth_inc_dexc", "dPhi_truth_inc_dexc", 100,-200,200);
  h_dPhi_truth_inc_dexc->GetXaxis()->SetTitle("dPhi_truth (deg)");
  
  h_dPhi_reco_inc = new TH1F("dPhi_reco_inc", "dPhi_reco_inc", 100,-200,200);
  h_dPhi_reco_inc->GetXaxis()->SetTitle("dPhi_reco (deg)");
  
  h_dPhi_truth_EvtMix_inc = new TH1F("dPhi_truth_EvtMix_inc", "dPhi_truth_EvtMix_inc", 100,-200,200);
  h_dPhi_truth_EvtMix_inc->GetXaxis()->SetTitle("dPhi_truth_EvtMix (deg)");
  
  h_dPhi_truth_corrected_inc = new TH1F("dPhi_truth_corrected_EvtMix_inc", "dPhi_truth_corrected_EvtMix_inc", 20,-200,200);
  h_dPhi_truth_corrected_inc->GetXaxis()->SetTitle("dPhi_truth (deg)");

  //------ for theta 70 to 90 degs
  h_dPhi_truth_70_90 = new TH1F("dPhi_truth_70_90", "dPhi_truth_70_90", 100,-200,200);
  h_dPhi_truth_70_90->GetXaxis()->SetTitle("dPhi_truth (deg)");

  h_dPhi_truth_70_90_dexc = new TH1F("dPhi_truth_70_90_dexc", "dPhi_truth_70_90_dexc", 100,-200,200);
  h_dPhi_truth_70_90_dexc->GetXaxis()->SetTitle("dPhi_truth (deg)");
  
  h_dPhi_reco_70_90 = new TH1F("dPhi_reco_70_90", "dPhi_reco_70_90", 100,-200,200);
  h_dPhi_reco_70_90->GetXaxis()->SetTitle("dPhi_reco (deg)");
  
  h_dPhi_truth_EvtMix_70_90 = new TH1F("dPhi_truth_EvtMix_70_90", "dPhi_truth_EvtMix_70_90", 100,-200,200);
  h_dPhi_truth_EvtMix_70_90->GetXaxis()->SetTitle("dPhi_truth_EvtMix (deg)");
  
  h_dPhi_truth_corrected_70_90 = new TH1F("dPhi_truth_corrected_EvtMix_70_90", "dPhi_truth_corrected_EvtMix_70_90", 20,-200,200);
  h_dPhi_truth_corrected_70_90->GetXaxis()->SetTitle("dPhi_truth (deg)");

  h_Edep_sum_reco = new TH1F("h_Edep_reco_sum", "h_Edep_reco_sum", 500,0,1.1);
  //h_test2d = new TH2F("test2D", "test2D", 80,-40,40, 80,-40,40);
  //h_test2d = new TH2F("test2D", "test2D", 10000,-200,200, 100,-200,200);
  h_test2d = new TH2F("test2D", "test2D", 100,0,1.1, 100,0,1.1);
}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* detType) {
  string nameDetType=detType;//vvv
  //TDirectory * dir = new TDirectory("TreeMaker2");
  TChain *tree = new TChain("tree");
  //  TChain *tree = new TChain("PreSelection");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of detType " << detType << std::endl;
  }

  if(nameDetType!="signalH") nameDetType="BG";
  if(nameDetType=="signalH") nameDetType="signal";
  //cout<<"Treating the input files as "<<nameDetType<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameDetType);
  
  BookHistogram(outFileName);

  //Jets = 0;
}
Bool_t AnalyzeLightBSM::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                   
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeLightBSM::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

AnalyzeLightBSM::~AnalyzeLightBSM() { 

  if (!fChain) {
    return;}
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

