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
 
//Defining the GAGG scintillator (Gd3Al2Ga3O12:Ce (Cerium is just activator which affects SY and decay timing. I am not adding it in the Material explicitly))
  
  G4Material *DetectorMat = new G4Material("LYSO", 7.1*g/cm3, 4);
  
  // G4Element *Gd = nist->FindOrBuildElement("G4_Gd");
  // G4Element *Al = nist->FindOrBuildElement("G4_Al");
  // G4Element *Ga = nist->FindOrBuildElement("G4_Ga");
  // G4Element *O = nist->FindOrBuildElement("G4_O");

  G4Element *Lu = nist->FindOrBuildElement("Lu");
  G4Element *Y = nist->FindOrBuildElement("Y");
  G4Element *Si = nist->FindOrBuildElement("Si");
  G4Element *O = nist->FindOrBuildElement("O");
  
  DetectorMat -> AddElement(Lu, 9);
  DetectorMat -> AddElement(Y, 1);
  DetectorMat -> AddElement(Si, 5);
  DetectorMat -> AddElement(O, 25);
  //G4cout << "\n\nName of the Detector Material is: " << DetectorMat->GetBaseMaterial()->GetName() << G4endl;
  
 

  std::vector<G4double> RIPhoEnergy = {0.996 *eV, 1.177 *eV, 1.425 *eV, 1.673 *eV, 1.968 *eV, 2.304 *eV, 2.643 *eV, 2.841 *eV, 2.949 *eV,
				       3.044 *eV, 3.130 *eV, 3.215 *eV, 3.328 *eV, 3.432 *eV, 3.680 *eV, 4.036 *eV, 4.383 *eV, 4.789 *eV};
  
  std::vector<G4double> RefractiveIdxDet = {1.784, 1.784, 1.789, 1.797, 1.805, 1.808, 1.816, 1.824, 1.833,
					    1.835, 1.833, 1.838, 1.835, 1.835, 1.841, 1.852, 1.868, 1.882};

  // std::vector<G4double> AbsPhoEnergy = {511 * keV};
  // std::vector<G4double> AbsLength = {10.4 * mm};

  std::vector<G4double> ScintPhoEnergy = {2.189 *eV, 2.199 *eV, 2.212 *eV, 2.222 *eV, 2.234 *eV, 2.254 *eV, 2.268 *eV, 2.285 *eV, 2.293 *eV, 2.316 *eV,
					  2.330 *eV, 2.350 *eV, 2.366 *eV, 2.373 *eV, 2.380 *eV, 2.396 *eV, 2.408 *eV, 2.411 *eV, 2.428 *eV, 2.444 *eV,
					  2.450 *eV, 2.458 *eV, 2.471 *eV, 2.493 *eV, 2.508 *eV, 2.524 *eV, 2.534 *eV, 2.538 *eV, 2.546 *eV, 2.573 *eV,
					  2.588 *eV, 2.593 *eV, 2.603 *eV, 2.619 *eV, 2.634 *eV, 2.650 *eV, 2.657 *eV, 2.668 *eV, 2.687 *eV, 2.714 *eV,
					  2.753 *eV, 2.784 *eV, 2.815 *eV, 2.824 *eV, 2.839 *eV, 2.880 *eV, 2.914 *eV, 2.924 *eV, 2.940 *eV, 2.940 *eV,
					  2.954 *eV, 2.962 *eV, 2.970 *eV, 2.984 *eV, 2.989 *eV, 3.001 *eV, 3.018 *eV, 3.028 *eV, 3.037 *eV, 3.090 *eV,
					  3.111 *eV, 3.144 *eV, 3.165 *eV, 3.195 *eV, 3.229 *eV, 3.247 *eV, 3.263 *eV, 3.285 *eV, 3.314 *eV, 3.354 *eV,
					  3.385 *eV, 3.400 *eV, 3.418 *eV, 3.435 *eV, 3.513 *eV, 3.559 *eV, 3.635 *eV};
  
  std::vector<G4double> ScintFastArray = {0.016, 0.004, 0.023, 0.002, 0.019, 0.026, 0.021, 0.019, 0.035, 0.033,
					  0.035, 0.047, 0.055, 0.039, 0.052, 0.055, 0.070, 0.064, 0.082, 0.084,
					  0.103, 0.100, 0.112, 0.127, 0.156, 0.167, 0.187, 0.204, 0.191, 0.244,
					  0.267, 0.284, 0.289, 0.323, 0.356, 0.365, 0.394, 0.411, 0.444, 0.504,
					  0.602, 0.698, 0.760, 0.791, 0.819, 0.887, 0.931, 0.941, 0.966, 0.951,
					  0.960, 0.990, 0.980, 0.990, 0.996, 0.980, 0.969, 0.934, 0.929, 0.803,
					  0.754, 0.658, 0.594, 0.501, 0.372, 0.309, 0.250, 0.185, 0.106, 0.050,
					  0.036, 0.008, 0.017, 0.006, 0.016, 0.006, 0.006};

  
  // std::vector<G4double> ScintPhoEnergy = {1.94 * eV, 2.51 * eV, 3.66 * eV};
  // std::vector<G4double> ScintFastArray = {0.019, 0.097, 0.020};
  // std::vector<G4double> ScintSlowArray = {0.171, 0.876, 0.183};

  // G4MaterialPropertyVector *MPVec;
  // MPVec->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintilSlowArray, false, true);

  //G4int lenArray = 10;
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
  
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, true);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 30000 / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40 * ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 1);
  myMPT1->AddConstProperty("SCINTILLATIONRISETIME1", 70 * ps);

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
  G4double rad_dz = 10 * cm;
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
 












































