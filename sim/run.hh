#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include <TROOT.h>
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

#include "G4Run.hh"
#include "detector.hh"

class MyRunAction : public G4UserRunAction
{
public:
  MyRunAction();
  ~MyRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  TFile *hfile;
  TTree *tree;
  //int x=90;
  std::vector<double> X, Y, Z, E;

  //tree->Branch("positionX", &p);
  
  //void FillTree(std::vector<G4double> m)
  //void FillTree()
  //void FillTree(G4double p)
  //{ //p = &m;
  //tree->Fill();
  //}

private:
  //int x=1;
  //std::vector<G4double> *Xcoord;
  
  
};

#endif

