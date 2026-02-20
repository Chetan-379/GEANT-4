//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  4 01:48:45 2016 by ROOT version 6.06/01
// from TTree PreSelection/PreSelection
// found on file: root://cmseos.fnal.gov//store/user/vhegde/myProduction_V11/Spring16.WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_96_RA2AnalysisTree.root
//////////////////////////////////////////////////////////
#ifndef NtupleVariables_h
#define NtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;
class NtupleVariables : public TSelector {
 public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
  Int_t           nOp_GAGG;
  Int_t           nOp_LYSO;
  Int_t           nOp_GAGG_QE;
  Int_t           nOp_LYSO_QE;
  Int_t           nOp_GAGG_truth;
  Int_t           nOp_LYSO_truth;
  Float_t         Edep_truth;
  Float_t         Edep_GAGG_truth;
  Float_t         Edep_LYSO_truth;
  Int_t           nOpGAGG_quad[4];
  Int_t           nOpLYSO_quad[4];
  
   // List of branches
  TBranch         *b_nOp_GAGG;   //!
  TBranch         *b_nOp_LYSO;   //!
  TBranch         *b_nOp_GAGG_QE;   //!
  TBranch         *b_nOp_LYSO_QE;   //!
  TBranch         *b_nOp_GAGG_truth;   //!
  TBranch         *b_nOp_LYSO_truth;   //!
  TBranch         *b_Edep_truth;   //!
  TBranch         *b_Edep_GAGG_truth;   //!
  TBranch         *b_Edep_LYSO_truth;   //!
  TBranch         *b_nOpGAGG_quad;   //!
  TBranch         *b_nOpLYSO_quad;   //!  
  
  
  
  NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   void    Init(TTree *tree, string);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  double  DeltaPhi(double, double);
  double  DeltaR(double eta1, double phi1, double eta2, double phi2);
  void    sortTLorVec(vector<TLorentzVector> *);
  double TransMass(double phi1, double phi2, double pt1, double pt2);
  double MinDr(TLorentzVector v1,vector<TLorentzVector> v2);
  double MinDr2(vector<TLorentzVector> v1,TLorentzVector v2);
  double getCrossSection(std::string process_name);
  std::map<std::string,float> cross_sectionValues;

};

#endif
#ifdef NtupleVariables_cxx
void NtupleVariables::Init(TTree *tree, string nameData)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

   // Set object pointer

   // left_module;
   
   // right_module[0][0] =0;

   // for (int i =0; i< 8; i++){
   //      for (int j =0; j< 8; j++){
   // 	  //for (int k =0; k< 8; k++){
   // 	  left_module[i][j] = Float_t(0.);
   // 	       right_module[i][j] = Float_t(0.);
   // 	       //}
   // 	}
   // }
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nOp_GAGG", &nOp_GAGG, &b_nOp_GAGG);
   fChain->SetBranchAddress("nOp_LYSO", &nOp_LYSO, &b_nOp_LYSO);
   fChain->SetBranchAddress("nOp_GAGG_QE", &nOp_GAGG_QE, &b_nOp_GAGG_QE);
   fChain->SetBranchAddress("nOp_LYSO_QE", &nOp_LYSO_QE, &b_nOp_LYSO_QE);
   fChain->SetBranchAddress("nOp_GAGG_truth", &nOp_GAGG_truth, &b_nOp_GAGG_truth);
   fChain->SetBranchAddress("nOp_LYSO_truth", &nOp_LYSO_truth, &b_nOp_LYSO_truth);
   fChain->SetBranchAddress("Edep_truth", &Edep_truth, &b_Edep_truth);
   fChain->SetBranchAddress("Edep_GAGG_truth", &Edep_GAGG_truth, &b_Edep_GAGG_truth);
   fChain->SetBranchAddress("Edep_LYSO_truth", &Edep_LYSO_truth, &b_Edep_LYSO_truth);
   fChain->SetBranchAddress("nOpGAGG_quad", &nOpGAGG_quad, &b_nOpGAGG_quad);
   fChain->SetBranchAddress("nOpLYSO_quad", &nOpLYSO_quad, &b_nOpLYSO_quad);
   Notify();
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef temp_cxx
