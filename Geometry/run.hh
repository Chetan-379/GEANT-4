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

  std::vector<G4double> HitEdep, HitTime, HitTrklen, HitScatAngle, HitEta, HitPol0, HitPol1, HitPol2, HitEout, HitEin, HitScatMomX, HitScatMomY, HitScatMomZ;
  std::vector<G4double> HitPosX, HitPosY, HitPosZ;
  std::vector<G4int> HitDetId, HitGunId, HitProcId;
 
private:
  
  
};

#endif
