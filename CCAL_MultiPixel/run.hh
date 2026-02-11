#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include <TROOT.h>
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "CellHit.hh"

#include "G4Run.hh"

class MyRunAction : public G4UserRunAction
{
public:
  MyRunAction();
  ~MyRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  TFile *hfile;
  TTree *tree;

  G4float Edep_truth, Edep_G_truth, Edep_L_truth;
  G4int nOpGAGG, nOpLYSO, nOpGAGG_truth, nOpLYSO_truth, nOpG_quad[4], nOpL_quad[4];

  // TH1F *h_Op_lmbda;
   TH2F *h_zVslmbda;
};

#endif
