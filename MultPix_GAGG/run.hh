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
  std::vector<double> X, Y, Z, E;
  G4double Total_E, Total_Compt_Edep, Photo_Edep;

  G4int nOpPhotons = 0;

  std::vector<G4double> HitEdep, HitTime, HitScatAngle, HitEta, HitEout, HitEin;
  std::vector<G4double> HitPosX, HitPosY, HitPosZ;
  std::vector<G4int> HitDetId, HitGunId, HitProcId;

  //std::vector<std::vector<G4int>> Module1, Module2;
  G4float Module1[8][8][6], Module2[8][8][6], Edep_truth, E_outPtcl=0., E_verify=0;

  G4int nPtcl_out =0, nSec_e=0, nSec_pho =0;

  std::vector<G4double> E_eOut, E_phoOut;

  G4float scat_theta;
  
  enum info_pack{nOpPho=0, DetPosX=1, DetPosY=2, nGenOp=3, Edep = 4, E_elec=5};  
  //CellHit Block1[8][8], Block2[8][8];
};

#endif
