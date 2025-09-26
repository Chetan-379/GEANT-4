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

  G4double nOpPhotons = 0, nOpPhotons_end =0;
  
private:
  MyRunAction *runObject;
  G4int ievent=1, chkEvt = 2693;
  G4int decade = 0;  
};
#endif
