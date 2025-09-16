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
  void FillHistogram(int, double, double, double, double, double, double, double, double, double, double, double);
  string MatName(int);
  void print(Long64_t);

  //define histograms here
  TH1D *h_Compt_Edep, *h_Photo_Edep, *h_Total_Edep, *h_Total_Edep_fine_binned;
  TH2D *h_ComptVsPhoto_Edep;
  TH2D *h_ScatAng_SDvsAna;
  TH1D *h_diff_Ana_SD;
  TH1I *h_nHits;
  TH1D *h_compSigEta[15], *h_hitTime[15];
  TH1D *h_theta[15], *h_eta[15];
  TH2D *h_theta_eta[15];
  TH1F *h_pol0[15], *h_pol1[15], *h_pol2[15];
  TH1F *h_hitPosX[15], *h_hitPosY[15], *h_hitPosZ[15];
  TH1F *h_Eout[15], *h_Ein[15];
  TH2D *h_theta_Eout[15];
  TH2F *h_posX_posY, *h_polX_polY;
  TH1F *h_Pol_Ana;
  TH2F *h_PolX_Ana_SD, *h_PolY_Ana_SD, *h_PolZ_Ana_SD;

  TH1I *h_nAnniPho;

  // TH1F *h_rough[5];
  // TH2F *h_rough2D[5];

  TH1F *h_dt_A1S1, *h_dt_A1S2, *h_dt_A2S1, *h_dt_A2S2, *h_dt_A1A2, *h_dt_S1S2;
  TH2F *h_STij;

  TH1F *h_dPhi_inVis,*h_dPhi_outVis, *h_dEta_inVis, *h_dEta_outVis, *h_dEta_inc;
  TH2F *h_Theta1vsTheta2, *h_Eta1VsEta2, *h_Eta1VsEta2_inc, *h_ThetaSimAna;
 
  TFile *oFile;
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //TGaxis::SetMaxDigits(2);

  h_Compt_Edep = new TH1D("Compton_Edep","Edep_via_Compton",100,0,1.1);
  h_Photo_Edep = new TH1D("Photo_Edep","Edep_via_Photoelecric",100,0,1.1);
   
  h_ComptVsPhoto_Edep = new TH2D("ComptVsPho","Compt_vs_Photo_Edep",100,0.,1.1, 100,0.,1.1);
  h_ComptVsPhoto_Edep->SetXTitle("Comp_Edep");
  h_ComptVsPhoto_Edep->SetYTitle("PhotoElectric_Edep");
    
  h_Total_Edep = new TH1D("Total_Edep","Total_Edep",100,0,1.1);
   
  h_ScatAng_SDvsAna = new TH2D("ScatAng_SDvsAna","ScatAng_SDvsAna",316,0.,3.16, 316,0.,3.16);
  h_ScatAng_SDvsAna->SetXTitle("Theta_Ana(rad)");
  h_ScatAng_SDvsAna->SetYTitle("Theta_SD(rad)");
   
  h_diff_Ana_SD = new TH1D("diff_Ana_SD", "diff_Ana_SD", 632,-3.16,3.16);
   
  h_posX_posY = new TH2F("XposVsYpos","XposVsYpos",1200,-600,600, 1200,-600,600);
  h_posX_posY->SetXTitle("X_Pos (mm)");
  h_posX_posY->SetYTitle("Y_Pos (mm)");
   
  h_polX_polY = new TH2F("XpolVsYpol","XpolVsYpol",120,-1.2,1.2, 120,-1.2,1.2);
  h_polX_polY->SetXTitle("X_Pol");
  h_polX_polY->SetYTitle("Y_Pol");
   
  h_PolX_Ana_SD = new TH2F("PolX_AnaVsSD","PolX_AnaVsSD",240,-1.2,1.2, 240,-1.2,1.2);
  h_PolX_Ana_SD->SetXTitle("PolX_Ana");
  h_PolX_Ana_SD->SetYTitle("PolX_SD");

  h_PolY_Ana_SD = new TH2F("PolY_AnaVsSD","PolY_AnaVsSD",240,-1.2,1.2, 240,-1.2,1.2);
  h_PolY_Ana_SD->SetXTitle("PolY_Ana");
  h_PolY_Ana_SD->SetYTitle("PolY_SD");
    
  h_PolZ_Ana_SD = new TH2F("PolZ_AnaVsSD","PolZ_AnaVsSD",240,-1.2,1.2, 240,-1.2,1.2);
  h_PolZ_Ana_SD->SetXTitle("PolZ_Ana");
  h_PolZ_Ana_SD->SetYTitle("PolZ_SD");
   
  h_Pol_Ana = new TH1F("Pol_Ana", "Pol_Ana", 240,-1.2,1.2);
   
  h_nAnniPho = new TH1I("nAnniPho", "nAnniPho", 10,0,10);

  h_dEta_inVis = new TH1F("dEta_inVis", "dEta_inVis", 400, 0, 400);
  h_dEta_outVis = new TH1F("dEta_outVis", "dEta_outVis", 400, 0, 400);
    
  h_dEta_inc = new TH1F("dEta_inc", "dEta_inc", 400, 0, 400);

  h_dPhi_inVis = new TH1F("DelPhi_highVis", "DelPhi_highVis", 100,-1,200);
  h_dPhi_outVis = new TH1F("DelPhi_lowVis", "DelPhi_lowVis", 100,-1,200);
    

  h_Eta1VsEta2 = new TH2F("Eta1VsEta2", "Eta1VsEta2", 100,-200,200, 100,-200,200);
  h_Eta1VsEta2->GetXaxis()->SetTitle("Eta1_(deg)");
  h_Eta1VsEta2->GetYaxis()->SetTitle("Eta2_(deg)");

  h_Eta1VsEta2_inc = new TH2F("Eta1VsEta2_inc", "Eta1VsEta2_inc", 100,-200,200, 100,-200,200);
  h_Eta1VsEta2_inc->GetXaxis()->SetTitle("Eta1_(deg)");
  h_Eta1VsEta2_inc->GetYaxis()->SetTitle("Eta2_(deg)");

  h_ThetaSimAna = new TH2F("Theta_SimAna", "Theta_SimAna", 100,0,200, 100,0,200);
  h_ThetaSimAna->GetXaxis()->SetTitle("thetaAna_(deg)");
  h_ThetaSimAna->GetYaxis()->SetTitle("thetaSim_(deg)");

  h_Theta1vsTheta2 = new TH2F("Theta1vsTheta2", "Theta1vsTheta2", 90,0,180, 90,0,180);
  h_Theta1vsTheta2->GetXaxis()->SetTitle("theta1_(deg)");
  h_Theta1vsTheta2->GetYaxis()->SetTitle("theta2_(deg)");
    
  h_dt_A1S1 = new TH1F("dt_A1S1", "dt_A1S1", 500, -5, 5);  //1000->500
  h_dt_A1S2 = new TH1F("dt_A1S2", "dt_A1S2", 500, -5, 5);
  h_dt_A2S1 = new TH1F("dt_A2S1", "dt_A2S1", 500, -5, 5);
  h_dt_A2S2 = new TH1F("dt_A2S2", "dt_A2S2", 500, -5, 5); 
  h_dt_A1A2 = new TH1F("dt_A1A2", "dt_A1A2", 500, -5, 5);
  h_dt_S1S2 = new TH1F("dt_S1S2", "dt_S1S2", 500, -5, 5);

  h_STij = new TH2F("h_ST1jVsST2j", "h_ST1jVsST2j", 100,-5,5, 100,-5,5);
  h_STij->GetXaxis()->SetTitle("ST1j (ns)");
  h_STij->GetYaxis()->SetTitle("ST2j (ns)");

  // h_rough[0] = new TH1F("rough0", "rough0", 4,0,4);
  // h_rough[0]->GetXaxis()->SetTitle("gId");

  // h_rough2D[0] = new TH2F("rough2D_0", "rough2D_0", 100,0,200, 100,0,200);
  // h_rough2D[0]->GetXaxis()->SetTitle("ThetaAna_(deg)");
  // h_rough2D[0]->GetYaxis()->SetTitle("ThetaSim_(deg)");


  double bin_edges[16];
  for (int j = 0; j <= 15; ++j) {
    bin_edges[j] = j - 0.5;
  }    
  h_nHits = new TH1I("nHits", "nHits", 15, bin_edges);  
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

