//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec 13 17:28:39 2025 by ROOT version 6.26/06
// from TTree tree/HitInfo
// found on file: src_root_files/test.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *positionX;
   vector<double>  *positionY;
   vector<double>  *positionZ;
   vector<double>  *Energy;
   Double_t        Total_Edep;
   Double_t        Edep_Compt;
   Double_t        Edep_Photo;
   Int_t           nOpticalPhotons;
   vector<double>  *Hit_Edep;
   vector<double>  *Hit_PositionX;
   vector<double>  *Hit_PositionY;
   vector<double>  *Hit_PositionZ;
   vector<double>  *Hit_Time;
   vector<int>     *Hit_DetId;
   vector<int>     *Hit_GunId;
   vector<int>     *Hit_ProcId;
   vector<double>  *Hit_ScatAngle;
   vector<double>  *Hit_Eta;
   vector<double>  *Hit_Eout;
   vector<double>  *Hit_Ein;
   Float_t         right_module[8][8][4];
   Float_t         left_module[8][8][4];

   // List of branches
   TBranch        *b_positionX;   //!
   TBranch        *b_positionY;   //!
   TBranch        *b_positionZ;   //!
   TBranch        *b_Energy;   //!
   TBranch        *b_Total_Edep;   //!
   TBranch        *b_Edep_Compt;   //!
   TBranch        *b_Edep_Photo;   //!
   TBranch        *b_nOpticalPhotons;   //!
   TBranch        *b_Hit_Edep;   //!
   TBranch        *b_Hit_PositionX;   //!
   TBranch        *b_Hit_PositionY;   //!
   TBranch        *b_Hit_PositionZ;   //!
   TBranch        *b_Hit_Time;   //!
   TBranch        *b_Hit_DetId;   //!
   TBranch        *b_Hit_GunId;   //!
   TBranch        *b_Hit_ProcId;   //!
   TBranch        *b_Hit_ScatAngle;   //!
   TBranch        *b_Hit_Eta;   //!
   TBranch        *b_Hit_Eout;   //!
   TBranch        *b_Hit_Ein;   //!
   TBranch        *b_Module1;   //!
   TBranch        *b_Module2;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("src_root_files/test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("src_root_files/test.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
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
   Hit_Edep = 0;
   Hit_PositionX = 0;
   Hit_PositionY = 0;
   Hit_PositionZ = 0;
   Hit_Time = 0;
   Hit_DetId = 0;
   Hit_GunId = 0;
   Hit_ProcId = 0;
   Hit_ScatAngle = 0;
   Hit_Eta = 0;
   Hit_Eout = 0;
   Hit_Ein = 0;
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
   fChain->SetBranchAddress("nOpticalPhotons", &nOpticalPhotons, &b_nOpticalPhotons);
   fChain->SetBranchAddress("Hit_Edep", &Hit_Edep, &b_Hit_Edep);
   fChain->SetBranchAddress("Hit_PositionX", &Hit_PositionX, &b_Hit_PositionX);
   fChain->SetBranchAddress("Hit_PositionY", &Hit_PositionY, &b_Hit_PositionY);
   fChain->SetBranchAddress("Hit_PositionZ", &Hit_PositionZ, &b_Hit_PositionZ);
   fChain->SetBranchAddress("Hit_Time", &Hit_Time, &b_Hit_Time);
   fChain->SetBranchAddress("Hit_DetId", &Hit_DetId, &b_Hit_DetId);
   fChain->SetBranchAddress("Hit_GunId", &Hit_GunId, &b_Hit_GunId);
   fChain->SetBranchAddress("Hit_ProcId", &Hit_ProcId, &b_Hit_ProcId);
   fChain->SetBranchAddress("Hit_ScatAngle", &Hit_ScatAngle, &b_Hit_ScatAngle);
   fChain->SetBranchAddress("Hit_Eta", &Hit_Eta, &b_Hit_Eta);
   fChain->SetBranchAddress("Hit_Eout", &Hit_Eout, &b_Hit_Eout);
   fChain->SetBranchAddress("Hit_Ein", &Hit_Ein, &b_Hit_Ein);
   fChain->SetBranchAddress("right_module", right_module, &b_Module1);
   fChain->SetBranchAddress("left_module", left_module, &b_Module2);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
