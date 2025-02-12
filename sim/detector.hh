#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
//#include "run.hh"
//#include "event.hh"


class MySensitiveDetector : public G4VSensitiveDetector
{
public:
  MySensitiveDetector(G4String);
  ~MySensitiveDetector();
  G4ThreeVector posDetector;
  G4double edep;
  

private:
  virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);
  //MyEventAction *fEventAction;
  
  
};

#endif

