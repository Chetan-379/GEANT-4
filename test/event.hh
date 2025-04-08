#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4AnalysisManager.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
public:
  MyEventAction(MyRunAction*);
  ~MyEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
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
  G4int n_Refraction, n_Reflection, n_TIR;
  
private:
  G4double fX;
  //std::vector<G4double> Xarray, Yarray, Zarray, Earray;
  MyRunAction *runObject;
  G4int ievent=1;
  G4double fEdep;
  G4double compt_total_edep, photo_total_edep;
  
  
};

#endif
