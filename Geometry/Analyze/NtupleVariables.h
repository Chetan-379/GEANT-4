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
   vector<double>  *positionX;
   vector<double>  *positionY;
   vector<double>  *positionZ;
   vector<double>  *Energy;
   Double_t        Total_Edep;
   Double_t        Edep_Compt;
   Double_t        Edep_Photo;
   // Int_t           nOpticalPhotons;
   // vector<double>  *OptPho_PosX;
   // vector<double>  *OptPho_PosY;
   // vector<double>  *OptPho_PosZ;
   // vector<double>  *OptPho_Energy;
   // vector<double>  *OptPho_Time;
  vector<double>  *Hit_Edep;
  vector<double>  *Hit_PositionX;
  vector<double>  *Hit_PositionY;
  vector<double>  *Hit_PositionZ;
  vector<double>  *Hit_Time;
  vector<double>  *Hit_TrkLen;
  vector<int>     *Hit_DetId;
  vector<string>  *Hit_ProcName;
  vector<int>  *Hit_ProcId;
  vector<double>  *Hit_ScatAngle;
  
   // List of branches
   TBranch        *b_positionX;   //!
   TBranch        *b_positionY;   //!
   TBranch        *b_positionZ;   //!
   TBranch        *b_Energy;   //!
   TBranch        *b_Total_Edep;   //!
   TBranch        *b_Edep_Compt;   //!
   TBranch        *b_Edep_Photo;   //!
   // TBranch        *b_nOpticalPhotons;
   // TBranch        *b_OptPho_PosX;
   // TBranch        *b_OptPho_PosY;
   // TBranch        *b_OptPho_PosZ;
   // TBranch        *b_OptPho_Energy;
   // TBranch        *b_OptPho_Time;
  TBranch         *b_Hit_Edep;
  TBranch         *b_Hit_PositionX;
  TBranch         *b_Hit_PositionY;
  TBranch         *b_Hit_PositionZ;
  TBranch         *b_Hit_Time;
  TBranch         *b_Hit_TrkLen;
  TBranch         *b_Hit_DetId;
  TBranch         *b_Hit_ProcName;
  TBranch         *b_Hit_ProcId;
  TBranch         *b_Hit_ScatAngle;
  
  
  
  NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   void    Init(TTree *tree, string);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  /*  double  DeltaPhi(double, double); */
  /*  double  DeltaR(double eta1, double phi1, double eta2, double phi2); */
  /*  void    sortTLorVec(vector<TLorentzVector> *);    */
  /*  double  MinDr(TLorentzVector, vector<TLorentzVector>); */
  /* double  MinDr2(vector<TLorentzVector>,TLorentzVector); */
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
   positionX = 0;
   positionY = 0;
   positionZ = 0;
   Energy = 0;
   // OptPho_PosX =0;
   // OptPho_PosY =0;
   // OptPho_PosZ =0;
   // OptPho_Energy =0;
   // OptPho_Time =0;
   Hit_Edep = 0;
   Hit_PositionX = 0;
   Hit_PositionY = 0;
   Hit_PositionZ = 0;
   Hit_Time = 0;
   Hit_TrkLen = 0;
   Hit_DetId = 0;
   //Hit_ProcName = 0;
   Hit_ProcName = 0;
   Hit_ScatAngle =0;
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("positionX", &positionX, &b_positionX);
   fChain->SetBranchAddress("positionY", &positionY, &b_positionY);
   fChain->SetBranchAddress("positionZ", &positionZ, &b_positionZ);
   fChain->SetBranchAddress("Energy", &Energy, &b_Energy);
   fChain->SetBranchAddress("Total_Edep", &Total_Edep, &b_Total_Edep);
   fChain->SetBranchAddress("Edep_Compt", &Edep_Compt, &b_Edep_Compt);
   fChain->SetBranchAddress("Edep_Photo", &Edep_Photo, &b_Edep_Photo);
   // fChain->SetBranchAddress("nOpticalPhotons", &nOpticalPhotons, &b_nOpticalPhotons);
   // fChain->SetBranchAddress("OptPho_PosX", &OptPho_PosX, &b_OptPho_PosX);
   // fChain->SetBranchAddress("OptPho_PosY", &OptPho_PosY, &b_OptPho_PosY);
   // fChain->SetBranchAddress("OptPho_PosZ", &OptPho_PosZ, &b_OptPho_PosZ);
   // fChain->SetBranchAddress("OptPho_Energy", &OptPho_Energy, &b_OptPho_Energy);
   // fChain->SetBranchAddress("OptPho_Time", &OptPho_Time, &b_OptPho_Time);
   fChain->SetBranchAddress("Hit_Edep", &Hit_Edep, &b_Hit_Edep);
   fChain->SetBranchAddress("Hit_PositionX", &Hit_PositionX, &b_Hit_PositionX);
   fChain->SetBranchAddress("Hit_PositionY", &Hit_PositionY, &b_Hit_PositionY);
   fChain->SetBranchAddress("Hit_PositionZ", &Hit_PositionZ, &b_Hit_PositionZ);
   fChain->SetBranchAddress("Hit_Time", &Hit_Time, &b_Hit_Time);
   fChain->SetBranchAddress("Hit_TrkLen", &Hit_TrkLen, &b_Hit_TrkLen);
   fChain->SetBranchAddress("Hit_DetId", &Hit_DetId, &b_Hit_DetId);
   fChain->SetBranchAddress("Hit_ProcName", &Hit_ProcName, &b_Hit_ProcName);
   fChain->SetBranchAddress("Hit_ProcId", &Hit_ProcId, &b_Hit_ProcId);
   fChain->SetBranchAddress("Hit_ScatAngle", &Hit_ScatAngle, &b_Hit_ScatAngle);

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
