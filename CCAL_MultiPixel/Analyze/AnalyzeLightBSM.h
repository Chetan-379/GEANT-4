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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"

#pragma link C++ class std::vector < std::vector> + ;
#pragma link C++ class std::vector < TLorentzVector> + ;

class AnalyzeLightBSM : public NtupleVariables
{

public:
  AnalyzeLightBSM(const TString &inputFileList = "foo.txt", const char *outFileName = "histo.root", const char *detType = "PbWO4");
  ~AnalyzeLightBSM();
  Bool_t FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void EventLoop(const char *, const char *, const char *);
  void BookHistogram(const char *);
  void FillHistogram(double, double, double, double, double, double, double, double);
  string MatName(int);
  void print(Long64_t);
  float GetScatAngle(float);
  float GetScatElecE(float);
  float Deg(float);
  float Rad(float);

  enum Quds{qd1 =0, qd2 =1, qd3 =2, qd4 =3};

  float Pi = TMath::Pi();
  
  // define histograms here 
  TH1I *h_nOp_GAGG, *h_nOp_LYSO, *h_nOp_GAGG_QE, *h_nOp_LYSO_QE, *h_nOp_GAGG_truth, *h_nOp_LYSO_truth, *h_nOp_GL_sum, *h_nOpG_q12, *h_nOpG_q13, *h_nOpL_q12, *h_nOpL_q13;

  TH1I  *h_nOpG_quad[4], *h_nOpL_quad[4];
  TH1F *h_Edep_truth, *h_EdepG_truth, *h_EdepL_truth;

  TH2F *h_nOpG_vs_Edep, *h_nOpL_vs_Edep, *h_nOpSumGL_vs_Edep, *h_nOpG_vs_EdepG, *h_nOpL_vs_EdepL, *h_nOpG_truth_vs_EdepG, *h_nOpL_truth_vs_EdepL, *h_nOpG_QE_vs_EdepG, *h_nOpL_QE_vs_EdepL;
  TH2I *h_nOp_GvsL, *h_nOp_GvsL_truth, *h_nOpG_q12_vs_q13, *h_nOpL_q12_vs_q13;

  TFile *oFile;
};

#endif

#ifdef AnalyzeLightBSM_cxx

void AnalyzeLightBSM::BookHistogram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);
  // TGaxis::SetMaxDigits(2);

  h_nOp_GAGG = new TH1I("n_Op_GAGG", "nOp_GAGG", 1500, 0 ,15000);
  h_nOp_GAGG ->GetXaxis()->SetTitle("nOp_GAGG");
  
  h_nOp_LYSO = new TH1I("n_Op_LYSO", "nOp_LYSO", 2000, 0 ,20000);
  h_nOp_LYSO ->GetXaxis()->SetTitle("nOp_LYSO");

  h_nOp_GAGG_QE = new TH1I("n_Op_GAGG_QE", "nOp_GAGG_QE", 1500, 0 ,15000);
  h_nOp_GAGG_QE ->GetXaxis()->SetTitle("nOp_GAGG_QE");
  
  h_nOp_LYSO_QE = new TH1I("n_Op_LYSO_QE", "nOp_LYSO_QE", 2000, 0 ,20000);
  h_nOp_LYSO_QE ->GetXaxis()->SetTitle("nOp_LYSO_QE");
  
  h_Edep_truth = new TH1F("Edep_truth", "Edep_truth", 100, 0 ,0.6);
  h_Edep_truth ->GetXaxis()->SetTitle("Edep_truth (MeV)");
  
  h_nOpG_vs_Edep = new TH2F("nOp_GAGG_vs_Edep", "nOp_GAGG_vs_Edep", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpG_vs_Edep->GetXaxis()->SetTitle("Edep (MeV)");
  h_nOpG_vs_Edep->GetYaxis()->SetTitle("nOp_rear_GAGG");

  h_nOpL_vs_Edep = new TH2F("nOp_LYSO_vs_Edep", "nOp_LYSO_vs_Edep", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpL_vs_Edep->GetXaxis()->SetTitle("Edep (MeV)");
  h_nOpL_vs_Edep->GetYaxis()->SetTitle("nOp_rear_LYSO");

  h_nOp_GvsL = new TH2I("nOp_GAGG_vs_LYSO", "nOp_GAGG_vs_LYSO", 1500, 0, 15000, 1500, 0, 15000);
  h_nOp_GvsL->GetXaxis()->SetTitle("nOp_GAGG");
  h_nOp_GvsL->GetYaxis()->SetTitle("nOp_LYSO");

  //true nOp-------
  h_nOp_GAGG_truth = new TH1I("n_Op_GAGG_truth", "nOp_GAGG_truth", 1500, 0 ,15000);
  h_nOp_GAGG_truth ->GetXaxis()->SetTitle("nOp_GAGG_truth");
  
  h_nOp_LYSO_truth = new TH1I("n_Op_LYSO_truth", "nOp_LYSO_truth", 1500, 0 ,15000);
  h_nOp_LYSO_truth ->GetXaxis()->SetTitle("nOp_LYSO_truth");

  h_nOp_GvsL_truth = new TH2I("nOp_GAGG_vs_LYSO_truth", "nOp_GAGG_vs_LYSO_truth", 1500, 0, 15000, 1500, 0, 15000);
  h_nOp_GvsL_truth->GetXaxis()->SetTitle("nOp_GAGG_truth");
  h_nOp_GvsL_truth->GetYaxis()->SetTitle("nOp_LYSO_truth");
  
  h_nOpG_truth_vs_EdepG = new TH2F("nOp_GAGG_truth_vs_EdepG", "nOp_GAGG_truth_vs_EdepG", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpG_truth_vs_EdepG->GetXaxis()->SetTitle("Edep_GAGG (MeV)");
  h_nOpG_truth_vs_EdepG->GetYaxis()->SetTitle("nOp_GAGG_truth");
  
  h_nOpL_truth_vs_EdepL = new TH2F("nOp_LYSO_truth_vs_EdepL", "nOp_LYSO_truth_vs_EdepL", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpL_truth_vs_EdepL->GetXaxis()->SetTitle("Edep_LYSO (MeV)");
  h_nOpL_truth_vs_EdepL->GetYaxis()->SetTitle("nOp_LYSO_truth"); 
  //--------------

  h_nOp_GL_sum = new TH1I("nOp_GAGG_LYSO_sum", "nOp_GAGG_LYSO_sum", 1500, 0 ,15000);
  h_nOp_GL_sum->GetXaxis()->SetTitle("nOp_sum(GAGG,LYSO)");

  h_nOpSumGL_vs_Edep = new TH2F("Sum_nOp_GL_vs_Edep_truth", "Sum_nOp_GL_vs_Edep_truth", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpSumGL_vs_Edep->GetXaxis()->SetTitle("Edep_truth (MeV)");
  h_nOpSumGL_vs_Edep->GetYaxis()->SetTitle("nOp_sum(GAGG,LYSO)");

  h_EdepG_truth = new TH1F("Edep_GAGG_truth", "Edep_GAGG_truth", 100, 0 ,0.6);
  h_EdepG_truth ->GetXaxis()->SetTitle("GAGG_Edep_truth (MeV)");

  h_EdepL_truth = new TH1F("Edep_LYSO_truth", "Edep_LYSO_truth", 100, 0 ,0.6);
  h_EdepL_truth ->GetXaxis()->SetTitle("LYSO_Edep_truth (MeV)");

  h_nOpG_vs_EdepG = new TH2F("nOp_GAGG_vs_Edep_GAGG", "nOp_GAGG_vs_Edep_GAGG", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpG_vs_EdepG->GetXaxis()->SetTitle("Edep_GAGG (MeV)");
  h_nOpG_vs_EdepG->GetYaxis()->SetTitle("nOp_GAGG");
  
  h_nOpL_vs_EdepL = new TH2F("nOp_LYSO_vs_Edep_LYSO", "nOp_LYSO_vs_Edep_LYSO", 100, 0, 0.6, 1500, 0, 15000);
  h_nOpL_vs_EdepL->GetXaxis()->SetTitle("Edep_LYSO (MeV)");
  h_nOpL_vs_EdepL->GetYaxis()->SetTitle("nOp_LYSO");

  h_nOpG_QE_vs_EdepG = new TH2F("nOp_GAGG_QE_vs_Edep_GAGG", "nOp_GAGG_QE_vs_Edep_GAGG", 100, 0, 0.6, 600, 0, 6000);
  h_nOpG_QE_vs_EdepG->GetXaxis()->SetTitle("Edep_GAGG (MeV)");
  h_nOpG_QE_vs_EdepG->GetYaxis()->SetTitle("nOp_GAGG_QE");
  
  h_nOpL_QE_vs_EdepL = new TH2F("nOp_LYSO_QE_vs_Edep_LYSO", "nOp_LYSO_QE_vs_Edep_LYSO", 100, 0, 0.6, 600, 0, 6000);
  h_nOpL_QE_vs_EdepL->GetXaxis()->SetTitle("Edep_LYSO (MeV)");
  h_nOpL_QE_vs_EdepL->GetYaxis()->SetTitle("nOp_LYSO_QE");



  //------------------
  vector<string> quads = {"quad1", "quad2", "quad3", "quad4"};
  char hname_nOpG_quad[50], hname_nOpL_quad[50];
  for (int i=0; i < quads.size(); i++){
    sprintf(hname_nOpG_quad, "%s_nOpGAGG",quads[i].c_str());
    sprintf(hname_nOpL_quad, "%s_nOpLYSO",quads[i].c_str());

    h_nOpG_quad[i] = new TH1I(hname_nOpG_quad, hname_nOpG_quad, 500,0,5000);
    h_nOpG_quad[i]->GetXaxis()->SetTitle(hname_nOpG_quad);
    
    h_nOpL_quad[i] = new TH1I(hname_nOpL_quad, hname_nOpL_quad, 500,0,5000);
    h_nOpL_quad[i]->GetXaxis()->SetTitle(hname_nOpL_quad);
  }

  h_nOpG_q12 = new TH1I ("nOp_GAGG_q12", "nOp_GAGG_q12", 1000,0,10000);
  h_nOpG_q12->GetXaxis()->SetTitle("nOp_GAGG_q12");

  h_nOpG_q13 = new TH1I ("nOp_GAGG_q13", "nOp_GAGG_q13", 1000,0,10000);
  h_nOpG_q13->GetXaxis()->SetTitle("nOp_GAGG_q13");
  
  h_nOpL_q12 = new TH1I ("nOp_LYSO_q12", "nOp_LYSO_q12", 1000,0,10000);
  h_nOpL_q12->GetXaxis()->SetTitle("nOp_LYSO_q12");
  
  h_nOpL_q13 = new TH1I ("nOp_LYSO_q13", "nOp_LYSO_q13", 1000,0,10000);
  h_nOpL_q13->GetXaxis()->SetTitle("nOp_LYSO_q13");
  
  h_nOpG_q12_vs_q13 = new TH2I("nOpG_q12_vs_q13", "nOpG_q12_vs_q13", 500,0,10000, 500,0,10000);
  h_nOpG_q12_vs_q13->GetXaxis()->SetTitle("nOp_GAGG_q12");
  h_nOpG_q12_vs_q13->GetYaxis()->SetTitle("nOp_GAGG_q13");

  h_nOpL_q12_vs_q13 = new TH2I("nOpL_q12_vs_q13", "nOpL_q12_vs_q13", 500,0,10000, 500,0,10000);
  h_nOpL_q12_vs_q13->GetXaxis()->SetTitle("nOp_LYSO_q12");
  h_nOpL_q12_vs_q13->GetYaxis()->SetTitle("nOp_LYSO_q13");
}

AnalyzeLightBSM::AnalyzeLightBSM(const TString &inputFileList, const char *outFileName, const char *detType)
{
  string nameDetType = detType; // vvv
  // TDirectory * dir = new TDirectory("TreeMaker2");
  TChain *tree = new TChain("tree");
  //  TChain *tree = new TChain("PreSelection");
  if (!FillChain(tree, inputFileList))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }
  else
  {
    std::cout << "Initiating analysis of detType " << detType << std::endl;
  }

  if (nameDetType != "signalH")
    nameDetType = "BG";
  if (nameDetType == "signalH")
    nameDetType = "signal";
  // cout<<"Treating the input files as "<<nameDetType<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree, nameDetType);

  BookHistogram(outFileName);

  // Jets = 0;
}
Bool_t AnalyzeLightBSM::FillChain(TChain *chain, const TString &inputFileList)
{

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if (!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while (1)
  {
    infile >> buffer;
    if (!infile.good())
      break;
    // std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeLightBSM::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (!fChain->InheritsFrom(TChain::Class()))
    return centry;
  TChain *chain = (TChain *)fChain;
  if (chain->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

AnalyzeLightBSM::~AnalyzeLightBSM()
{

  if (!fChain)
  {
    return;
  }
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}

#endif
