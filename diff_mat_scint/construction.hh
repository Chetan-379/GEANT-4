#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>
#include "G4LogicalBorderSurface.hh"
#include "G4GenericMessenger.hh"

#include "detector.hh"
#include "globals.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }
  
  virtual G4VPhysicalVolume *Construct();
  MySensitiveDetector *sensDet;

private:  
  virtual void ConstructSDandField();

  G4LogicalVolume *fScoringVolume;

  G4GenericMessenger *fMessenger;
  G4double depth;

  G4LogicalVolume *logicWorld, *logicDetector, *logicEncVol;

  G4Box *solidWorld, *solidEncVol;
  G4Trd *solidDetector;

  G4VPhysicalVolume *physWorld, *physDetector, *physEncVol;

  G4Material *worldMat, *DetectorMat;
  G4Element *Gd, *Al, *Ga, *O, *Lu, *Y, *Si, *La, *Br;

  void DefineMaterials();
};
  
#endif
