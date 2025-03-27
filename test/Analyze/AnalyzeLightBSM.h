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

#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;

class AnalyzeLightBSM : public NtupleVariables{

 public:
  AnalyzeLightBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *detType="PbWO4");
  //  std::cout<<"alpana"<<std::endl;
  ~AnalyzeLightBSM();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *,const char *);
  void     BookHistogram(const char *);
  void print(Long64_t);

  //define histograms here
  //TH1F *h_selectBaselineYields_;
  TH1D *h_Compt_Edep, *h_Photo_Edep, *h_Total_Edep, *h_Total_Edep_fine_binned;
  TH2D *h_ComptVsPhoto_Edep;
  TH1I *h_nOptPho, *h_OptPhoOnDet;
  TH2D *h_nOptPho_Edep;
  TH2I *h_nOptPhoOnDet_genOptPho;
  TH1D *h_OptPho_lmbda;

  TFile *oFile;
  /* TH1F *h_events; */
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //Initialize histogram here
  //h_selectBaselineYields_ = new TH1F("cutflows","cutflows",60,-0.5,60.5);
  h_Compt_Edep = new TH1D("Compton_Edep","Edep_via_Compton",70,0,0.7);
  h_Photo_Edep = new TH1D("Photo_Edep","Edep_via_Photoelecric",70,0,0.7);

  h_ComptVsPhoto_Edep = new TH2D("ComptVsPho","Compt_vs_Photo_Edep",70,0.,0.7, 70,0.,0.7);
  h_ComptVsPhoto_Edep->SetXTitle("Compton_Edep");
  h_ComptVsPhoto_Edep->SetYTitle("PhotoElectric_Edep");

  h_Total_Edep = new TH1D("Total_Edep","Total_Edep",70,0,0.7);
  h_Total_Edep_fine_binned = new TH1D("Edep_fine", "Edep_fine", 700, 0.45, 0.52);

  h_nOptPho = new TH1I("nOptical_Photons", "nOptical_Photons", 60, 0, 60);

  h_nOptPho_Edep = new TH2D("nOptPhoVsEdep", "nOptPhoVsEdep", 70, 0., 0.7, 60, 0, 60);
  h_nOptPho_Edep->SetXTitle("Total Edep");
  h_nOptPho_Edep->SetYTitle("nOptical Photons");

  h_OptPhoOnDet = new TH1I("nOptPhoton_OnDet", "nOptPhoton_OnDet", 50, 0, 50);

  h_nOptPhoOnDet_genOptPho = new TH2I("OptPhotOnDet_vs_genOptPhot", "OptPhotOnDet_vs_genOptPhot", 50, 0,50, 50, 0, 50);
  h_nOptPhoOnDet_genOptPho -> SetXTitle("nOptPhotOnDet");
  h_nOptPhoOnDet_genOptPho -> SetYTitle("genOptPhotons");

  h_OptPho_lmbda = new TH1D("OptPho_lmbda", "OptPho_lmbda",500,0,1000);

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

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

