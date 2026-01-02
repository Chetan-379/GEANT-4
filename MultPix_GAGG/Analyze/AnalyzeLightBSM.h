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

  TH1I *h_nOpPho_crys, *h_nOpPho_gen_crys, *h_nOpPho_mdl, *h_nOpPho_gen_mdl;
  TH1F *h_DetX_crys, *h_Edep_crys, *h_Edep_mdl;

  TH2D *h_nOp_vs_Edep_crys, *h_nOp_vs_Edep_mdl;

  TH1F* h_Edep_diff, *h_Edep_scat, *h_Edep_abs;

  TH2F* h_Edep_true_Vs_reco_crys, *h_Edep_scatVsabs;

  TH1I* h_nHits_mdl;

  TH2F* h_DetXvsDetY;

  TH1F* h_theta_reco, *h_theta_truth;

  TH2F* h_elecE_vs_Edep_crys;
   
  TFile *oFile;

  //enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep=4};
  enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep = 4, E_elec=5};
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

  h_DetX_crys = new TH1F("h_DetPosX_crys", "h_DetPosX_crys",1000,-40,40);

  h_nOpPho_crys = new TH1I("h_nOpPho_cryst", "h_nOpPho_cryst",3000,0,30000);
  h_nOpPho_gen_crys = new TH1I("h_nOpPho_gen_cryst", "h_nOpPho_gen_cryst",3000,0,30000);
  
  h_nOpPho_mdl = new TH1I("h_nOpPho_mdl", "h_nOpPho_mdl",3000,0,30000);
  h_nOpPho_gen_mdl = new TH1I("h_nOpPho_gen_mdl", "h_nOpPho_gen_mdl",3000,0,30000);
  

  h_Edep_crys = new TH1F ("h_Edep_cryst", "h_Edep_cryst", 100,0,1.1);
  h_nOp_vs_Edep_crys = new TH2D("nOp_vs_Edep_cryst", "nOp_vs_Edep_cryst", 100,0.,1.1, 3000,0,30000);

  h_Edep_mdl = new TH1F ("h_Edep_mdl", "h_Edep_mdl", 100,0,1.1);
  h_nOp_vs_Edep_mdl = new TH2D("nOp_vs_Edep_mdl", "nOp_vs_Edep_mdl", 100,0.,1.1, 3000,0,30000);

  h_Edep_diff = new TH1F ("h_Edep_diff", "h_Edep_diff", 1000,0.,0.5);

  h_Edep_true_Vs_reco_crys = new TH2F("Edep_TrueVsReco_cryst", "Edep_TrueVsReco_cryst", 100,0.,1.1, 100,0,1.1);
  h_Edep_true_Vs_reco_crys->GetXaxis()->SetTitle("Edep_{truth} (MeV)");
  h_Edep_true_Vs_reco_crys->GetYaxis()->SetTitle("Edep_{reco} (MeV)");

  h_nHits_mdl = new TH1I("h_nHits_mdl", "h_nHits_mdl", 10,0,10);

  h_DetXvsDetY = new TH2F("DetXvsDetY", "DetXvsDetY", 100,-40,40, 100,-40,40);

  h_theta_reco = new TH1F("h_theta_reco", "h_theta_reco", 400,0,200);
  h_theta_truth = new TH1F("h_theta_truth", "h_theta_truth", 400,0,200);

  h_elecE_vs_Edep_crys= new TH2F("scat_ElecVsEdep", "scat_ElecVsEdep", 100,0.,1.1, 100,0.,1.1);
  h_elecE_vs_Edep_crys->GetXaxis()->SetTitle("Edep_truth (MeV)");
  h_elecE_vs_Edep_crys->GetYaxis()->SetTitle("E_scat_elec (MeV)");

  h_Edep_scat = new TH1F("Edep_scat_truth", "Edep_scat_truth", 100,0,0.8);
  h_Edep_abs = new TH1F("Edep_abs_truth", "Edep_abs_truth", 100,0,0.8);

  h_Edep_scatVsabs = new TH2F("Edep_scatVsabs", "Edep_scatVsabs", 100,0,0.8, 100,0,0.8);
  h_Edep_scatVsabs->GetXaxis()->SetTitle("scat_Edep (MeV)");
  h_Edep_scatVsabs->GetYaxis()->SetTitle("abs_Edep (MeV)");
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

