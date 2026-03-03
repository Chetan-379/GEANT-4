#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpBoundaryProcess.hh"


#include "construction.hh"
#include "event.hh"

class MySteppingAction : public G4UserSteppingAction
{
public:
  MySteppingAction(MyEventAction* eventAction, MyRunAction* runAction);
  ~MySteppingAction();

  virtual void UserSteppingAction(const G4Step*);

private:
  MyEventAction *fEventAction;
  MyRunAction *fRunAction;
  
  G4OpBoundaryProcess* fBoundary = nullptr;  
};

#endif
