//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 24 12:11:26 2025 by ROOT version 6.28/10
// from TTree tree/HitInfo
// found on file: compt_photo_PbWO4_10000evts.root
//////////////////////////////////////////////////////////

#ifndef tree_h
#define tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class tree {
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
   Double_t        diff_Edep_ComptPhoto;

   // List of branches
   TBranch        *b_positionX;   //!
   TBranch        *b_positionY;   //!
   TBranch        *b_positionZ;   //!
   TBranch        *b_Energy;   //!
   TBranch        *b_Total_Edep;   //!
   TBranch        *b_Edep_Compt;   //!
   TBranch        *b_Edep_Photo;   //!
   TBranch        *b_diff_Edep_ComptPhoto;   //!

   tree(TTree *tree=0);
   virtual ~tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree_cxx
tree::tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("compt_photo_PbWO4_10000evts.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("compt_photo_PbWO4_10000evts.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
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

void tree::Init(TTree *tree)
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
   fChain->SetBranchAddress("diff_Edep_ComptPhoto", &diff_Edep_ComptPhoto, &b_diff_Edep_ComptPhoto);
   Notify();
}

Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx
