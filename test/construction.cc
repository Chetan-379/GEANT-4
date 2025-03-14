#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  G4NistManager *nist = G4NistManager::Instance();

  G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexAerogel[2] = {1.1, 1.1};
  G4double rindexWorld[2] = {1.0, 1.0};
   
  G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);

  G4Isotope *Na22 = new G4Isotope("Na22", 11, 22, 21.9944 * g / mole);
  G4Element *elNa22 = new G4Element("Sodium-22", "Na22", 1);
  elNa22-> AddIsotope(Na22, 100.0 * perCent);
  
  G4Material *matNa22 = new G4Material("Na22Source", 0.97 * g / cm3, 1);
  matNa22->AddElement(elNa22, 100.0 * perCent);

  G4Box *solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");

  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld, "physWorld", 0, false, 0, true);

  //spherical Radioactive source
  G4double sourceRadius = 1. * mm;
  
  G4Sphere *solidSource = new G4Sphere("solidSource", 0.0, sourceRadius, 0.0, 360. * deg, 0.0, 180. * deg);
  G4LogicalVolume *logicSource = new G4LogicalVolume(solidSource, matNa22, "logicSource");
  G4VPhysicalVolume * physSource = new G4PVPlacement(0, G4ThreeVector(0., 0., -1. * cm), logicSource, "physSource", logicWorld, false, 0, true);

  G4VisAttributes *sourceVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 1.0, 0.5));
  sourceVisAtt->SetForceSolid(true);
  logicSource->SetVisAttributes(sourceVisAtt);
  
  // Trapezoid shape                                                                                     
  G4double rad_dxa = 20 * cm, rad_dxb = 20 * cm;  
  G4double rad_dz = 5.0 * cm;
  G4double rad_dya = (rad_dxa+2*rad_dz), rad_dyb = (rad_dxa + 2 * rad_dz);
  
  G4Trd *solidDetector =
    new G4Trd("solidDetector",  // its name
  	      0.5 * rad_dxa, 0.5 * rad_dya, 0.5 * rad_dxb, 0.5 * rad_dyb,
              0.5 * rad_dz);  // its size	      

  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_PbWO4");
  G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_BGO");
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_Si");

   // G4Element *Cd = nist->FindOrBuildElement("G4_Cd");
   // G4Element *Zn = nist->FindOrBuildElement("G4_Zn");
   // G4Element *Te = nist->FindOrBuildElement("G4_Te");

  //  G4Material *DetectorMat = new G4Material("CdZnTe", 5.76*g/cm3, 3);
  // G4cout << "\n\nworking till here\n\n" << G4endl;
  //  DetectorMat->AddElement(Cd, 40*perCent);
  //   G4cout << "\n\nworking till here\n\n" << G4endl;
  //  DetectorMat->AddElement(Zn, 10*perCent);
  //   G4cout << "\n\nworking till here\n\n" << G4endl;
  //  DetectorMat->AddElement(Te, 50*perCent);
  //   G4cout << "\n\nworking till here\n\n" << G4endl;

   

  DetectorMat->SetMaterialPropertiesTable(mptWorld);

  logicDetector = new G4LogicalVolume(solidDetector,  // its solid
  						      DetectorMat, // its material
  						      "logicalDetector");  // its name

  fScoringVolume = logicDetector;
  
  G4RotationMatrix new_rotmX, new_rotmY, new_rotmZ;
     
  G4int Pos[2] = {-1,1};

  new_rotmX.rotateX(90 * deg);
  new_rotmY.rotateY(90 * deg);  
  new_rotmZ.rotateZ(90 * deg);

  G4ThreeVector position;
  
  for (G4int i=0; i<3; i++){     
    for (G4int j=0; j<2; j++){
      std::vector<G4RotationMatrix> new_rotm = {new_rotmX, new_rotmY, new_rotmZ};

      if (i==0) position = G4ThreeVector(0.0, Pos[j]*((rad_dxa+rad_dz)/2), 0.0);
      if (i==1) position = G4ThreeVector(-1*Pos[j]*((rad_dxa+rad_dz)/2), 0.0, 0.0);      
      if (i==2) position = G4ThreeVector(0.0, 0.0, -1*Pos[j]*((rad_dxa+rad_dz)/2));
      
     
      G4Transform3D transform = G4Transform3D(new_rotm[i], position);
      if (i==0) new_rotmX.rotateX(180 * deg);
      if (i==1) new_rotmY.rotateY(180 * deg);    
      if (i==2) new_rotmZ.rotateX(180 * deg);
    
      G4VPhysicalVolume *physDetector = new G4PVPlacement(transform, logicDetector, "physDetector", logicWorld, false, (i+1)*(j+1), true);      
    }
  }
 

  //============================Working for cuboidal volumes Don't Edit===================================================================================================
  // G4Box *solidRadiator = new G4Box("solidRadiator", 0.1*m, 0.1*m, 0.025*m);

  // G4Material *RadiatorMat = nist->FindOrBuildMaterial("G4_Si");

  // RadiatorMat->SetMaterialPropertiesTable(mptWorld);

  
  // G4LogicalVolume *logicRadiator = new G4LogicalVolume(solidRadiator, RadiatorMat, "logicalRadiator");


  
  // G4RotationMatrix new_rotmX, new_rotmY, new_rotmZ;
     
  // G4int Pos[2] = {-1,1};

  // new_rotmX.rotateX(90 * deg);
  // new_rotmY.rotateY(90 * deg);  
  // new_rotmZ.rotateZ(90 * deg);

  // std::vector<G4RotationMatrix> new_rotm = {new_rotmX, new_rotmY, new_rotmZ};
  // G4ThreeVector position;
  
  // for (G4int i=0; i<3; i++){     
  //   for (G4int j=0; j<2; j++){
      
  //     if (i==0) position = G4ThreeVector(0.0,Pos[j]*0.125*m, 0.0);
  //     if (i==1) position = G4ThreeVector(Pos[j]*0.125*m,0.0, 0.0);
  //     if (i ==2) position = G4ThreeVector(0.0, 0.0, Pos[j]*0.125*m);
     
  //     G4Transform3D transform = G4Transform3D(new_rotm[i], position);
  //     G4VPhysicalVolume *physRadiator = new G4PVPlacement(transform, logicRadiator, "physRadiator", logicWorld, false, (i+1)*(j+1), true);
  //}
  //}
  
  //===============================================================================================================================================
     
   
    return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
  sensDet = new MySensitiveDetector("SensitiveDetector");
  logicDetector->SetSensitiveDetector(sensDet);
}
