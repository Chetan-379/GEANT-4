#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
//#include <TROOT.h>
#include "G4AnalysisManager.hh"

#include "run.hh"

class MyEventAction : public G4UserEventAction
{
public:
  MyEventAction(MyRunAction*);
  ~MyEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  //void AddX(G4double XCoord) {fX += XCoord;}
  void AddX(G4double XCoord) {fX += 1;}
  // void StoreHitPos(G4ThreeVector Position)
  // {   
  //   HitsArray.push_back(Position);
  // }

  void StoreX(G4double XCoord)
  {
    Xarray.push_back(XCoord);
  }

  //int x=90;
  //MyRunAction runObject;
  // void StoreZ(G4double ZCoord)
  // {
  //   Zarray.push_back(ZCoord);
  // }
  
  //int x=1;  

private:
  G4double fX;
  int x=56;
  std::vector<G4double> Xarray; //, Yarray, Zarray;
  //std::vector<G4ThreeVector> HitsArray;
  MyRunAction *runObject;
  

// public:
//   void FillTree()
//   {
//     //
//     //runObject.tree->Fill();
//   }

};

#endif
