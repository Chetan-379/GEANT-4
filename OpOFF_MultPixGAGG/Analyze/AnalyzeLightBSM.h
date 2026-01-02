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

  TH1I *h_nOpPho_crys[2], *h_nOpPho_gen_crys[2], *h_nOpPho_mdl[2], *h_nOpPho_gen_mdl[2];
  TH1F *h_DetX_crys[2], *h_Edep_crys[2], *h_Edep_mdl[2];

  TH2D *h_nOp_vs_Edep_crys[2], *h_nOp_vs_Edep_mdl[2];

  TH1F* h_Edep_diff[2], *h_Edep_scat[2], *h_Edep_abs[2];

  TH2F* h_Edep_true_Vs_reco_crys[2], *h_Edep_scatVsabs[2];

  TH1I* h_nHits_mdl[2];

  TH2F* h_DetXvsDetY[2];

  TH1F* h_theta_reco[2], *h_theta_truth[2], *h_scatElec_E[2], *h_phi_truth[2], *h_dPix[2];

  TH2F* h_elecE_vs_Edep_crys[2], *h_elecE_vs_scatEdep[2];

  TH1F *h_dPhi_truth, *h_dPhi_rmBkg_truth;

  TH2F *h_Phi1_vs_Phi2_truth;
   
  TFile *oFile;

  //enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep=4};
  enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep = 4, E_elec=5, tryVar =6};
};
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
  
  char hname_DetPosX_crys[50], hname_nOpPho_cryst[50], hname_nOpPho_gen_cryst[50], hname_nOpPho_mdl[50], hname_nOpPho_gen_mdl[50], hname_Edep_cryst[50], hname_nOp_vs_Edep_cryst[50], hname_Edep_mdl[50], hname_nOp_vs_Edep_mdl[50], hname_Edep_diff[50], hname_Edep_TrueVsReco_cryst[50], hname_nHits_mdl[50], hname_DetXvsDetY[50], hname_theta_reco[50], hname_theta_truth[50], hname_scat_ElecVsEdep[50], hname_scat_ElecVSscat_Edep[50], hname_Edep_scat_truth[50], hname_Edep_abs_truth[50], hname_Edep_scatVsabs[50], hname_scat_elec_E[50], hname_phi_truth[50], hname_dPix[50];

  for (int i=0; i<2; i++){
  sprintf(hname_DetPosX_crys, "%s_DetPosX_crys", ModName[i].c_str());
  sprintf(hname_nOpPho_cryst, "%s_nOpPho_cryst",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_cryst, "%s_nOpPho_gen_cryst",ModName[i].c_str());
  sprintf(hname_nOpPho_mdl, "%s_nOpPho_mdl",ModName[i].c_str());
  sprintf(hname_nOpPho_gen_mdl, "%s_nOpPho_gen_mdl",ModName[i].c_str());
  sprintf(hname_Edep_cryst, "%s_Edep_cryst",ModName[i].c_str());
  sprintf(hname_nOp_vs_Edep_cryst, "%s_nOp_vs_Edep_cryst",ModName[i].c_str());
  sprintf(hname_Edep_mdl, "%s_Edep_mdl",ModName[i].c_str());
  sprintf(hname_nOp_vs_Edep_mdl, "%s_nOp_vs_Edep_mdl",ModName[i].c_str());
  sprintf(hname_Edep_diff, "%s_Edep_diff",ModName[i].c_str());
  sprintf(hname_Edep_TrueVsReco_cryst, "%s_Edep_TrueVsReco_cryst",ModName[i].c_str());
  sprintf(hname_nHits_mdl, "%s_nHits_mdl",ModName[i].c_str());
  sprintf(hname_DetXvsDetY, "%s_DetXvsDetY",ModName[i].c_str());
  sprintf(hname_theta_reco, "%s_theta_reco",ModName[i].c_str());
  sprintf(hname_theta_truth, "%s_theta_truth",ModName[i].c_str());
  sprintf(hname_scat_ElecVsEdep, "%s_scat_ElecVsEdep",ModName[i].c_str());
  sprintf(hname_scat_ElecVSscat_Edep, "%s_scat_ElecVSscat_Edep",ModName[i].c_str());
  sprintf(hname_Edep_scat_truth, "%s_Edep_scat_truth",ModName[i].c_str());
  sprintf(hname_Edep_abs_truth, "%s_Edep_abs_truth",ModName[i].c_str());
  sprintf(hname_Edep_scatVsabs, "%s_Edep_scatVsabs",ModName[i].c_str());
  sprintf(hname_scat_elec_E, "%s_scat_elec_E",ModName[i].c_str());
  sprintf(hname_phi_truth, "%s_phi_truth",ModName[i].c_str());
  sprintf(hname_dPix, "%s_pixel_distance",ModName[i].c_str());
  

  
  h_DetX_crys[i] = new TH1F(hname_DetPosX_crys, hname_DetPosX_crys,1000,-40,40);

  h_nOpPho_crys[i] = new TH1I(hname_nOpPho_cryst, hname_nOpPho_cryst,3000,0,30000);
  h_nOpPho_gen_crys[i] = new TH1I(hname_nOpPho_gen_cryst, hname_nOpPho_gen_cryst,3000,0,30000);
  
  h_nOpPho_mdl[i] = new TH1I(hname_nOpPho_mdl, hname_nOpPho_mdl,3000,0,30000);
  h_nOpPho_gen_mdl[i] = new TH1I(hname_nOpPho_gen_mdl, hname_nOpPho_gen_mdl,3000,0,30000);
  

  h_Edep_crys[i] = new TH1F (hname_Edep_cryst, hname_Edep_cryst, 100,0,1.1);
  h_nOp_vs_Edep_crys[i] = new TH2D(hname_nOp_vs_Edep_cryst, hname_nOp_vs_Edep_cryst, 100,0.,1.1, 3000,0,30000);

  h_Edep_mdl[i] = new TH1F (hname_Edep_mdl, hname_Edep_mdl, 100,0,1.1);
  h_nOp_vs_Edep_mdl[i] = new TH2D(hname_nOp_vs_Edep_mdl, hname_nOp_vs_Edep_mdl, 100,0.,1.1, 3000,0,30000);

  h_Edep_diff[i] = new TH1F (hname_Edep_diff, hname_Edep_diff, 1000,0.,0.5);

  h_Edep_true_Vs_reco_crys[i] = new TH2F(hname_Edep_TrueVsReco_cryst, hname_Edep_TrueVsReco_cryst, 100,0.,1.1, 100,0,1.1);
  h_Edep_true_Vs_reco_crys[i]->GetXaxis()->SetTitle("Edep_{truth} (MeV)");
  h_Edep_true_Vs_reco_crys[i]->GetYaxis()->SetTitle("Edep_{reco} (MeV)");

  h_nHits_mdl[i] = new TH1I(hname_nHits_mdl, hname_nHits_mdl, 10,0,10);

  h_DetXvsDetY[i] = new TH2F(hname_DetXvsDetY, hname_DetXvsDetY, 100,-40,40, 100,-40,40);

  h_theta_reco[i] = new TH1F(hname_theta_reco, hname_theta_reco, 400,0,200);
  h_theta_truth[i] = new TH1F(hname_theta_truth, hname_theta_truth, 400,0,200);

  h_elecE_vs_Edep_crys[i] = new TH2F(hname_scat_ElecVsEdep, hname_scat_ElecVsEdep, 100,0.,1.1, 100,0.,1.1);
  h_elecE_vs_Edep_crys[i]->GetXaxis()->SetTitle("Edep_truth (MeV)");
  h_elecE_vs_Edep_crys[i]->GetYaxis()->SetTitle("E_scat_elec (MeV)");

  h_elecE_vs_scatEdep[i] = new TH2F(hname_scat_ElecVSscat_Edep, hname_scat_ElecVSscat_Edep, 100,0.,1.1, 100,0.,1.1);
  h_elecE_vs_scatEdep[i]->GetXaxis()->SetTitle("Edep_scat (MeV)");
  h_elecE_vs_scatEdep[i]->GetYaxis()->SetTitle("E_scat_elec (MeV)");


  h_Edep_scat[i] = new TH1F(hname_Edep_scat_truth, hname_Edep_scat_truth, 100,0,0.8);
  h_Edep_abs[i] = new TH1F(hname_Edep_abs_truth, hname_Edep_abs_truth, 100,0,0.8);

  h_Edep_scatVsabs[i] = new TH2F(hname_Edep_scatVsabs, hname_Edep_scatVsabs, 100,0,0.8, 100,0,0.8);
  h_Edep_scatVsabs[i]->GetXaxis()->SetTitle("scat_Edep (MeV)");
  h_Edep_scatVsabs[i]->GetYaxis()->SetTitle("abs_Edep (MeV)");

  h_scatElec_E[i]=new TH1F(hname_scat_elec_E, hname_scat_elec_E, 100,0,0.8);

  h_phi_truth[i] = new TH1F(hname_phi_truth, hname_phi_truth, 100,-200,200);

  h_dPix[i] = new TH1F(hname_dPix, hname_dPix, 40,0,40);
  }
  //---------------------------
  
  
  h_dPhi_truth = new TH1F("diff_Phi_truth", "diff_Phi_truth", 100,-200,200);
  h_dPhi_truth->GetXaxis()->SetTitle("dPhi (deg)");
  
  h_Phi1_vs_Phi2_truth = new TH2F("Phi1_vs_Phi2_truth", "Phi1_vs_Phi2_truth", 100,-200,200, 100,-200,200);
  h_Phi1_vs_Phi2_truth->GetXaxis()->SetTitle("Phi1 (deg)");
  h_Phi1_vs_Phi2_truth->GetYaxis()->SetTitle("Phi2 (deg)");
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
  cout<<"Treating the input files as "<<nameDetType<<" for setting tree branches"<<endl;
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

