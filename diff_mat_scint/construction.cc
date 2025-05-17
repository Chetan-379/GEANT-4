                //=================================construction file for plastic scintillator============================ 

#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  //-----------------Defining and Adding the Material Properties----------
  //
  
  G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexWorld[2] = {1.0, 1.0};

  //World(Air)
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);

  //Detector Tile
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_PbWO4");
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_BGO");
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_Si");

  //Defining plastic Scintillator Material EJ230 which is PolyVinylToluene (C9H10)
  G4Material *DetectorMat = new G4Material("C9H10", 1.032*g/cm3, 2);
  G4Element *C = nist->FindOrBuildElement("C");
  G4Element *H = nist->FindOrBuildElement("H");

  DetectorMat->AddElement(C, 91.5*perCent);
  DetectorMat->AddElement(H, 8.5*perCent);
  
  // std::vector<G4double> RIPhoEnergy = {2.58 * eV, 2.95 * eV, 4.00 * eV};
  // std::vector<G4double> RefractiveIdxDet = {2.15, 2.2, 2.39};

  std::vector<G4double> RIPhoEnergy = {2.58 * eV, 2.95 * eV, 4.00 * eV};
  std::vector<G4double> RefractiveIdxDet = {1.58, 1.58, 1.58};


  std::vector<G4double> AbsPhoEnergy = {511 * keV};
  std::vector<G4double> AbsLength = {10.4 * mm};
  
  // std::vector<G4double> ScintPhoEnergy={1.95 * eV, 2.01 * eV, 2.08 * eV, 2.24 * eV, 2.51 * eV, 2.59 * eV, 2.69 * eV, 3.09 * eV, 3.40 * eV, 3.66 * eV};
  // std::vector<G4double> ScintFastArray={0.0213, 0.038, 0.0517, 0.0727, 0.0973, 0.0913, 0.078, 0.05, 0.0333, 0.02};
  // std::vector<G4double> ScintSlowArray={0.19, 0.34, 0.47, 0.65, 0.88, 0.82, 0.70, 0.45, 0.30, 0.18};


  std::vector<G4double> ScintPhoEnergy = {2.51 * eV, 2.52 * eV, 2.53 * eV, 2.55 * eV, 2.56 * eV, 2.57 * eV, 2.59 * eV, 2.60 * eV, 2.62 * eV, 2.63 * eV,
					  2.65 * eV, 2.66 * eV, 2.68 * eV, 2.69 * eV, 2.71 * eV, 2.72 * eV, 2.74 * eV, 2.76 * eV, 2.77 * eV, 2.79 * eV,
					  2.80 * eV, 2.82 * eV, 2.84 * eV, 2.86 * eV, 2.87 * eV, 2.89 * eV, 2.91 * eV, 2.93 * eV, 2.94 * eV, 2.96 * eV,
					  2.98 * eV, 3.00 * eV, 3.02 * eV, 3.04 * eV, 3.06 * eV, 3.08 * eV, 3.10 * eV, 3.12 * eV, 3.14 * eV, 3.16 * eV,
					  3.18 * eV, 3.20 * eV, 3.23 * eV, 3.25 * eV, 3.27 * eV, 3.29 * eV, 3.32 * eV, 3.34 * eV, 3.37 * eV, 3.39 * eV,
					  3.41 * eV, 3.44 * eV, 3.46 * eV, 3.49 * eV, 3.52 * eV};
    

  std::vector<G4double> ScintFastArray = {0.008, 0.015, 0.015, 0.015, 0.023, 0.023, 0.031, 0.031, 0.038, 0.046,
					  0.046, 0.054, 0.062, 0.069, 0.077, 0.085, 0.100, 0.108, 0.123, 0.138,
					  0.162, 0.185, 0.208, 0.238, 0.269, 0.300, 0.331, 0.369, 0.415, 0.462,
					  0.515, 0.577, 0.646, 0.708, 0.785, 0.862, 0.915, 0.962, 0.985, 1.000,
					  1.000, 0.992, 0.962, 0.908, 0.785, 0.638, 0.477, 0.308, 0.200, 0.131,
					  0.077, 0.046, 0.031, 0.015, 0.008};


  

  // std::vector<G4double> ScintPhoEnergy = {1.94 * eV, 2.51 * eV, 3.66 * eV};
  // std::vector<G4double> ScintFastArray = {0.019, 0.097, 0.020};
  // std::vector<G4double> ScintSlowArray = {0.171, 0.876, 0.183};

  // G4MaterialPropertyVector *MPVec;
  // MPVec->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintilSlowArray, false, true);

  G4int lenArray = 10;
  auto myMPT1 = new G4MaterialPropertiesTable();
  
  myMPT1->AddProperty("RINDEX", RIPhoEnergy, RefractiveIdxDet);

  // Check that group velocity is calculated from RINDEX
  if (myMPT1->GetProperty("RINDEX")->GetVectorLength()
      != myMPT1->GetProperty("GROUPVEL")->GetVectorLength())
  {
    G4ExceptionDescription ed;
    ed << "Error calculating group velocities. Incorrect number of entries "
          "in group velocity material property vector.";
    G4Exception("OpNovice::OpNoviceDetectorConstruction", "OpNovice001", FatalException, ed);
  }

    // Adding a property from two std::vectors. Argument createNewKey is false
  // and spline is true.

    
  //myMPT1->AddProperty("ABSLENGTH", AbsPhoEnergy, AbsLength, false, true);


  //G4cout << "\n\ndatatype is: " <<typeid( myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, true)).name() << G4endl;

  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, true);
  //myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", ScintPhoEnergy, ScintSlowArray, false, true);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 9700. / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1.5 * ns);
  // myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 300. * ns); 
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 1);
  //myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.9);
  myMPT1->AddConstProperty("SCINTILLATIONRISETIME1", 0.5 * ns);

  // Set the Birks Constant for scintillator
  DetectorMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
  
  //==============================================
  DetectorMat->SetMaterialPropertiesTable(myMPT1);
  //==============================================
    

  //----------Volumes-------------
  //
  // The world
  G4Box *solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");
  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld, "physWorld", 0, false, 0, true);

  // The enclosure
  G4Box *solidEncVol = new G4Box("EncVol", 10*cm, 10*cm, 10*cm);
  G4LogicalVolume *logicEncVol = new G4LogicalVolume(solidEncVol, worldMat, "logicalEncVol");
  G4VPhysicalVolume *physEncVol = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicEncVol, "physEncVol", logicWorld, false, 0, true);

  
  //Trapezoidal enclosure tiles                                                          
  G4double rad_dxa = 20 * cm, rad_dxb = 20 * cm;  
  G4double rad_dz = 20 * cm;
  G4double rad_dya = (rad_dxa+2*rad_dz), rad_dyb = (rad_dxa + 2 * rad_dz);

  G4RotationMatrix new_rotmX, new_rotmY, new_rotmZ;
  G4int Pos[2] = {-1,1};
  new_rotmX.rotateX(90 * deg);
  new_rotmY.rotateY(90 * deg);  
  new_rotmZ.rotateZ(90 * deg);
  
  G4ThreeVector position;
  
  G4Trd *solidDetector = new G4Trd("solidDetector",  // its name
				   0.5 * rad_dxa, 0.5 * rad_dya, 0.5 * rad_dxb, 0.5 * rad_dyb,
				   0.5 * rad_dz);  // its size	      

  logicDetector = new G4LogicalVolume(solidDetector,  // its solid
				      DetectorMat, // its material
				      "logicalDetector");  // its name

  // G4cout << "\n\n=====================working====================================\n\n" << G4endl;
  //Loops to place the tiles
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
  
  fScoringVolume = logicDetector;
  return physWorld;
  
}

void MyDetectorConstruction::ConstructSDandField()
{
  sensDet = new MySensitiveDetector("SensitiveDetector");
  logicDetector->SetSensitiveDetector(sensDet);  
}




//==================BACK UP===================================

//-----------MAking CZT Material---------------
//G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_BGO");
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


  // ============================ Working for adding the grid (Don't edit or remove)===================================================================================================
  // G4Box *solidDetector = new G4Box("solidDetector", 2.5* cm, 2.5* cm, 2.5* cm);
  // G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_Si");

  // DetectorMat->SetMaterialPropertiesTable(mptWorld);
  
  // logicDetector = new G4LogicalVolume(solidDetector, DetectorMat, "logicalDetector");
  
  // //G4RotationMatrix new_rotmX, new_rotmY, new_rotmZ;
     
  // G4int Pos[2] = {-1,1};

  // // new_rotmX.rotateX(90 * deg);
  // // new_rotmY.rotateY(90 * deg);  
  // // new_rotmZ.rotateZ(90 * deg);

  // //std::vector<G4RotationMatrix> new_rotm = {new_rotmX, new_rotmY, new_rotmZ};
  // G4ThreeVector position;
  
  // for (G4int s=0; s<3; s++){     
  //   for (G4int i=0; i<2; i++){
  //     for (G4int j=0; j<4; j++){
  // 	for (G4int k=0; k<4; k++){
  // 	  if(s==0) position = G4ThreeVector((-10 + (2.5*((2*j)+1))) * cm, Pos[i]*12.5 * cm, (-10 + (2.5*((2*k)+1))) * cm);
  // 	  if(s==1) position = G4ThreeVector(Pos[i]*12.5 * cm, (-10 + (2.5*((2*j)+1))) * cm, (-10 + (2.5*((2*k)+1))) * cm);
  // 	  if(s==2) position = G4ThreeVector((-10 + (2.5*((2*j)+1))) * cm, (-10 + (2.5*((2*k)+1))) * cm, Pos[i]*12.5 * cm);
  // 	  // if (i==1) position = G4ThreeVector(Pos[j]*0.125*m,0.0, 0.0);
  // 	  // if (i ==2) position = G4ThreeVector(0.0, 0.0, Pos[j]*0.125*m);
	
  // 	  //G4Transform3D transform = G4Transform3D(new_rotm[i], position);
  // 	  G4VPhysicalVolume *physDetector = new G4PVPlacement(0, position, logicDetector, "physDetector", logicWorld, false, (j+1)*(k+1), true);
  // 	}
  //     }
  //   }
  // }
  
  

 // ============================Working for cuboidal volumes Don't Edit===================================================================================================
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
  //   }
    
  // }
  
  //===============================================================================================================================================


//-------------------Definining Radioactive Source Volume------------
 // //spherical Radioactive source
  // G4double sourceRadius = 1. * mm;
  
  // G4Sphere *solidSource = new G4Sphere("solidSource", 0.0, sourceRadius, 0.0, 360. * deg, 0.0, 180. * deg);
  // G4LogicalVolume *logicSource = new G4LogicalVolume(solidSource, matNa22, "logicSource");
  // G4VPhysicalVolume * physSource = new G4PVPlacement(0, G4ThreeVector(0., 0., 0. * cm), logicSource, "physSource", logicWorld, false, 0, true);

  // G4VisAttributes *sourceVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 1.0, 0.5));
  // sourceVisAtt->SetForceSolid(true);
  // logicSource->SetVisAttributes(sourceVisAtt);
 












































