#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include <TROOT.h>
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"

#include "G4Run.hh"

class MyRunAction : public G4UserRunAction
{
public:
  MyRunAction();
  ~MyRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  G4double OpPho_lmbda, OpPho_CreateZ;

  TTree *tree;
 
private:
  TFile *hfile;
  
};

#endif
