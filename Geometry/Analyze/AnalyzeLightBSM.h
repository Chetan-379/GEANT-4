#ifndef AnalyzeLightBSM_H
#define AnalyzeLightBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
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
  //AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *detType="PbWO4");
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *detType="PbWO4");
  //  std::cout<<"alpana"<<std::endl;
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *);
  void     BookHistogram(const char *);
  void print(Long64_t);

  //define histograms here
  TH1D *h_Compt_Edep, *h_Photo_Edep, *h_Total_Edep, *h_Total_Edep_fine_binned;
  TH2D *h_ComptVsPhoto_Edep;
  TH2D *h_ScatAng_SDvsAna;
  TH1D *h_diff_Ana_SD;
  TH1I *h_nHits;
  TH1D *h_compSigEta;
  TH1D *h_theta, *h_eta;
  TH2D *h_theta_eta;
  TH1F *h_pol0;
  TH1F *h_pol1;
  TH1F *h_pol2;

  TFile *oFile;
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //TGaxis::SetMaxDigits(2);
  //Initialize histogram here
  //h_selectBaselineYields_ = new TH1F("cutflows","cutflows",60,-0.5,60.5);

  // vector<TString> Energy, Material;
  // vector<double> hist_Max;
  // vector<int> Energy_val;
  
  // Material = {"BGO", "PbWO4", "Plastic"};
  // Energy = {"511keV", "100keV", "150keV", "300keV", "450keV", "600keV", "800keV", "1000keV"};
  // Energy_val = {511, 100, 150, 300, 450, 600, 800, 1000};
  
  
  h_Compt_Edep = new TH1D("Compton_Edep","Edep_via_Compton",100,0,1.1);
  h_Photo_Edep = new TH1D("Photo_Edep","Edep_via_Photoelecric",100,0,1.1);

  h_ComptVsPhoto_Edep = new TH2D("ComptVsPho","Compt_vs_Photo_Edep",100,0.,1.1, 100,0.,1.1);
  h_ComptVsPhoto_Edep->SetXTitle("Comp_Edep");
  h_ComptVsPhoto_Edep->SetYTitle("PhotoElectric_Edep");

  h_Total_Edep = new TH1D("Total_Edep","Total_Edep",100,0,1.1);

  h_ScatAng_SDvsAna = new TH2D("ScatAng_SDvsAna","ScatAng_SDvsAna",316,0.,3.16, 316,0.,3.16);
  h_ScatAng_SDvsAna->SetXTitle("Theta_Ana(rad)");
  h_ScatAng_SDvsAna->SetYTitle("Theta_SD(rad)");

  h_nHits = new TH1I("nHits", "nHits", 15, 0, 15);
  h_diff_Ana_SD = new TH1D("diff_Ana_SD","diff_Ana_SD",632,-3.16,3.16);

  h_theta = new TH1D("scat_theta", "scat_theta", 200,-200, 200);
  h_theta->GetXaxis()->SetTitle("theta (deg)");
  
  h_eta = new TH1D("pol_scat_eta", "pol_scat_eta", 200, -200, 200);
  h_eta->GetXaxis()->SetTitle("eta (deg)");

  h_compSigEta = new TH1D("Eta_compSig", "Eta_compSig",200,-200,200);
  h_compSigEta->GetXaxis()->SetTitle("eta (deg)");

  h_theta_eta = new TH2D("ThetaVsEta", "ThetaVsEta", 200, -200, 200, 200, -200, 200);
  h_theta_eta->GetXaxis()->SetTitle("eta (deg)");
  h_theta_eta->GetYaxis()->SetTitle("theta (deg)");

  h_pol0 = new TH1F("Pol_0", "Pol_0", 2000,-1.2,1.2);
  h_pol1 = new TH1F("Pol_1", "Pol_1", 2000,-1.2,1.2);
  h_pol2 = new TH1F("Pol_2", "Pol_2", 2000,-1.2,1.2);
  
 



  //h_Total_Edep_fine_binned = new TH1D("Edep_fine", "Edep_fine", 700, 0.45, 0.52);

  //h_nOptPho = new TH1I("nOptical_Photons", "nOptical_Photons", 35000, 0, 35000);
  //h_nOptPho = new TH1I("nOptical_Photons", "nOptical_Photons", 67000, 0, 67000);
  //h_nOptPho = new TH1I("nOptical_Photons", "nOptical_Photons", 250, 0, 250);

  // h_nOptPho_Edep = new TH2D("nOptPhoVsEdep", "nOptPhoVsEdep", 70, 0., 0.7, 60, 0, 60);
  // h_nOptPho_Edep->SetXTitle("Total Edep");
  // h_nOptPho_Edep->SetYTitle("nOptical Photons");

  // h_OptPhoOnDet = new TH1I("nOptPhoton_OnDet", "nOptPhoton_OnDet", 50, 0, 50);

  // h_nOptPhoOnDet_genOptPho = new TH2I("OptPhotOnDet_vs_genOptPhot", "OptPhotOnDet_vs_genOptPhot", 50, 0,50, 50, 0, 50);
  // h_nOptPhoOnDet_genOptPho -> SetXTitle("nOptPhotOnDet");
  // h_nOptPhoOnDet_genOptPho -> SetYTitle("genOptPhotons");

  //h_OptPho_lmbda = new TH1D("OptPho_lmbda", "OptPho_lmbda",2000,0,2000);
  
  //h_OptPho_time = new TH1D("OptPho_time", "OptPho_time",5000,0,1000);
  //h_OptPho_time = new TH1D("OptPho_time", "OptPho_time",5000,0,1000);

  // h_OptPho_PosX = new TH1D("OptPho_PosX", "OptPho_PosX",400,-200,200);
  // h_OptPho_PosY = new TH1D("OptPho_PosY", "OptPho_PosY",400,-200,200);
  // h_OptPho_PosZ = new TH1D("OptPho_PosZ", "OptPho_PosZ",600,-300,300);
  
  // h_OptPho_XvsY = new TH2F("OptPho_XvsY", "OptPho_XvsY",400, -200, 200, 400, -200, 200);
  // h_OptPho_XvsY->SetXTitle("PosX");
  // h_OptPho_XvsY->SetYTitle("PosY");

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

