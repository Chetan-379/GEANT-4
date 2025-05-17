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

  TFile *hfile;
  TTree *tree;
  std::vector<double> X, Y, Z, E;
  G4double Total_E, Total_Compt_Edep, Photo_Edep;
  G4int nOpticalPhotons, nOptPhoOnDetEnd;
  std::vector<G4double> Opt_Photon_PosX, Opt_Photon_PosY, Opt_Photon_PosZ, Opt_Photon_Energy, Opt_Photon_Time;
  
 
private:
  
  
};

#endif

