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
  //TH3F *h_posXYZ;
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

    // h_posXYZ = new TH3F("PosXYZ","PosXYZ",1200,-600,600, 1200,-600,600, 1000, -500, 500);
    // h_posXYZ->SetXTitle("X (mm)");
    // h_posXYZ->SetYTitle("Y (mm)");
    // h_posXYZ->SetZTitle("Z (mm)");
    

     double bin_edges[16];
    for (int j = 0; j <= 15; ++j) {
      bin_edges[j] = j - 0.5;
    }    
    h_nHits = new TH1I("nHits", "nHits", 15, bin_edges);



    //vector<string> nHit_cat = {"all_inc", "nHit_eq1", "nHit_gt1", "nHit_eq2", "nHit_gt2"};
    vector<string> nHit_cat = {"all_inc", "nHit0", "nHit1_S", "nHit1_B", "nHit2_S1S1", "nHit2_S1S2", "nHit2_SB", "nHit2_B1B1", "nHit2_B1B2", "nHit2_BS", "remain", "test"};

  char hname_Compt_Edep[100], hname_Photo_Edep[100], hname_Total_Edep[100], hname_Total_Edep_fine_binned[100], hname_ComptVsPhoto_Edep[10], hname_ScatAng_SDvsAna[100], hname_diff_Ana_SD[100], hname_nHits[100], hname_compSigEta[100], hname_hitTime[100], hname_theta[100], hname_eta[100], hname_theta_eta[100], hname_pol0[100], hname_pol1[100], hname_pol2[100], hname_hitPosX[100], hname_hitPosY[100], hname_hitPosZ[100], hname_Eout[100], hname_Ein[100], hname_theta_Eout[100];

  for(int i=0; i< nHit_cat.size(); i++){            
    // sprintf(hname_h_Compt_Edep, "h_Compt_Edep_%s",nHit_cat[i].c_str());
    // sprintf(hname_Photo_Edep, "h_ComptVsPhoto_Edep_%s",nHit_cat[i].c_str());
    // sprintf(hname_Total_Edep, "h_Total_Edep_%s",nHit_cat[i].c_str());
    // sprintf(hname_Total_Edep_fine_binned, "h_Total_Edep_fine_binned_%s",nHit_cat[i].c_str());
    // sprintf(hname_ComptVsPhoto_Edep, "h_ComptVsPhoto_Edep_%s",nHit_cat[i].c_str());
    // sprintf(hname_ScatAng_SDvsAna, "h_ScatAng_SDvsAna_%s",nHit_cat[i].c_str());
    // sprintf(hname_diff_Ana_SD, "h_diff_Ana_SD_%s",nHit_cat[i].c_str());
    
    sprintf(hname_nHits, "h_nHits_%s",nHit_cat[i].c_str());
    sprintf(hname_compSigEta, "h_compSigEta_%s",nHit_cat[i].c_str());
    sprintf(hname_hitTime, "h_hitTime_%s",nHit_cat[i].c_str());
    sprintf(hname_theta, "h_theta_%s",nHit_cat[i].c_str());
    sprintf(hname_eta, "h_eta_%s",nHit_cat[i].c_str());
    sprintf(hname_theta_eta, "h_theta_eta_%s",nHit_cat[i].c_str());
    sprintf(hname_pol0, "h_pol0_%s",nHit_cat[i].c_str());
    sprintf(hname_pol1, "h_pol1_%s",nHit_cat[i].c_str());
    sprintf(hname_pol2, "h_pol2_%s",nHit_cat[i].c_str());
    sprintf(hname_hitPosX, "h_hitPosX_%s",nHit_cat[i].c_str());
    sprintf(hname_hitPosY, "h_hitPosY_%s",nHit_cat[i].c_str());
    sprintf(hname_hitPosZ, "h_hitPosZ_%s",nHit_cat[i].c_str());
    sprintf(hname_Eout, "h_Eout_%s",nHit_cat[i].c_str());
    sprintf(hname_Ein, "h_Ein_%s",nHit_cat[i].c_str());
    sprintf(hname_theta_Eout, "h_theta_Eout_%s", nHit_cat[i].c_str());
    
       
   
    h_hitTime[i] = new TH1D(hname_hitTime, hname_hitTime, 100,0 ,10);
  
    h_hitPosX[i] = new TH1F(hname_hitPosX, hname_hitPosX, 250,-700 ,700);
    h_hitPosY[i] = new TH1F(hname_hitPosY, hname_hitPosY, 250,-700 ,700);
    h_hitPosZ[i] = new TH1F(hname_hitPosZ, hname_hitPosZ, 250,-700 ,700);      
    
    h_theta[i] = new TH1D(hname_theta, hname_theta, 100,0, 200);
    h_theta[i]->GetXaxis()->SetTitle("theta (deg)");
  
    h_eta[i] = new TH1D(hname_eta, hname_eta, 200, -200, 200);
    h_eta[i]->GetXaxis()->SetTitle("eta (deg)");
    
    h_compSigEta[i] = new TH1D(hname_compSigEta, hname_compSigEta, 200,-200,200);
    h_compSigEta[i]->GetXaxis()->SetTitle("eta (deg)");
    
    h_theta_eta[i] = new TH2D(hname_theta_eta, hname_theta_eta, 200, -200, 200, 100, 0, 200);

    h_theta_eta[i]->GetXaxis()->SetTitle("eta (deg)");
    h_theta_eta[i]->GetYaxis()->SetTitle("theta (deg)");
    
    h_pol0[i] = new TH1F(hname_pol0, hname_pol0, 2000,-1.2,1.2);
    h_pol1[i] = new TH1F(hname_pol1, hname_pol1, 2000,-1.2,1.2);
    h_pol2[i] = new TH1F(hname_pol2, hname_pol2, 2000,-1.2,1.2);
    
    h_Eout[i] = new TH1F(hname_Eout, hname_Eout, 260, 0, 0.520);
    h_Eout[i]->GetXaxis()->SetTitle("Eout (MeV)");

    h_Ein[i] = new TH1F(hname_Ein, hname_Ein, 260, 0, 0.520);
    h_Ein[i]->GetXaxis()->SetTitle("Ein (MeV)");

    h_theta_Eout[i] = new TH2D(hname_theta_Eout, hname_theta_Eout, 500,0,200, 260,0.,0.520);
    h_theta_Eout[i]->GetYaxis()->SetTitle("Eout (MeV)");
    h_theta_Eout[i]->GetXaxis()->SetTitle("theta (deg)");        
  }
  
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

