#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4AnalysisManager.hh"
#include "G4AttValue.hh"

#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"

#include "run.hh"
#include "CellHit.hh"
#include "globals.hh"


class MyEventAction : public G4UserEventAction
{
public:
  MyEventAction(MyRunAction*);
  ~MyEventAction();

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);
  bool Compt_Before_Photo = false;
  G4int nCompt_Before_Photo;
  std::vector<G4double> Compt_edep, Photo_edep;  
  void AddX(G4double XCoord) {fX += 1;}

  void StorePos(G4double XCoord, G4double YCoord, G4double ZCoord, G4double edep)
  {
    Xarray.push_back(XCoord);
    Yarray.push_back(YCoord);
    Zarray.push_back(ZCoord);
    Earray.push_back(edep);
  }

  void AddEdep(G4double edep) {fEdep += edep;}

  std::vector<G4double> Xarray, Yarray, Zarray, Earray;
  std::vector<G4double> OptPho_PosX, OptPho_PosY, OptPho_PosZ, OptPho_Energy, OptPho_time;
  G4int nOptPho, TrkOnDet;
  
private:
  G4double fX;
  MyRunAction *runObject;
  G4int ievent=1;
  G4double fEdep;
  G4double compt_total_edep, photo_total_edep;

  CellHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
  void PrintEventStatistics(G4double CellEdep, G4double CellTrackLength) const;

  G4int fInStripHCID = -1;
  G4int fStripHCID = -1;

  std::vector<G4double> HitEdep_vec, HitTime_vec, HitTrkLen_vec;
  std::vector<double> HitPosX_vec, HitPosY_vec, HitPosZ_vec;
  std::vector<G4int> HitDetId_vec;
};

#endif
