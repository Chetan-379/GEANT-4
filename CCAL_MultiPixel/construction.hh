#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PVReplica.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "globals.hh"




class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  MyDetectorConstruction();
  ~MyDetectorConstruction();
  
  virtual G4VPhysicalVolume *Construct();

private:
  virtual void ConstructSDandField();
  G4Material *worldMat, *gagg, *lyso, *SiPM_Mat;;

  G4Element *Gd, *Al, *Ga, *O, *Lu, *Y, *Si;

  G4OpticalSurface* OpSurface_ScintAir, *OpSurface_SiPM;

  void DefineMaterials();
  
 
};
  
#endif
