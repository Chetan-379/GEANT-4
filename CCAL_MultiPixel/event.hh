#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4AnalysisManager.hh"
#include "G4AttValue.hh"

#include "G4RunManager.hh"
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

  G4int ievent=0;
  G4int chkEvt = 13010;

  std::vector<G4bool> store_log;
  std::vector<G4int> Det_row_Idx, Det_col_Idx, Det_clr_Idx;

  G4bool storeHit;
  
private:
  MyRunAction *runObject;
};

#endif
