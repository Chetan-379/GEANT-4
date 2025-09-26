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
  string MatName(int);
  void print(Long64_t);

  TH1F *h_CreateZ, *h_OpPho_lmbda;
  TH2F *h_lmbdVsZ;

  TFile *oFile;
};
#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName) {
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  //TGaxis::SetMaxDigits(2);

  h_CreateZ = new TH1F("OpPho_CreationZ", "OpPho_CreationZ", 100,-10,10);
  h_OpPho_lmbda = new TH1F("OpPho_lmbda", "OpPho_lmbda", 300,200,800);
  h_lmbdVsZ = new TH2F("OpPho_LamdaVsCreationZ", "OpPho_LamdaVsCreationZ", 200,-10,10, 300,200,800);
  h_lmbdVsZ->GetXaxis()->SetTitle("CreatePosZ(mm)");
  h_lmbdVsZ->GetYaxis()->SetTitle("OpPho_lambda(nm)");
}


AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char* detType) {
  string nameDetType=detType;//vvv
  //TDirectory * dir = new TDirectory("TreeMaker2");
  //TChain *tree = new TChain("OpPho_Ana");
  TChain *tree = new TChain("OpPho_End");
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

