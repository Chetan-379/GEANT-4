#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "G4AnalysisManager.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
public:
  MyEventAction(MyRunAction*);
  ~MyEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void AddX(G4double XCoord) {fX += XCoord;}
  // vector<G4double> StoreX(G4double XCoord)
  // {
  //   std::vector<G4double> Xarray;
  //   Xarray.push_back(fX);
  //   return Xarray;
  // }

private:
  G4double fX;
  
};

#endif
