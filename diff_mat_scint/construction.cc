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
  G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_BGO");
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_Si");
  
  
  // Material properties can be added as arrays. However, in this case it is
  // up to the user to make sure both arrays have the same number of elements.

  std::vector<G4double> RIPhoEnergy = {2.58 * eV, 2.95 * eV, 4.00 * eV};
  std::vector<G4double> RefractiveIdxDet = {2.15, 2.2, 2.39};

  std::vector<G4double> AbsPhoEnergy = {511 * keV};
  std::vector<G4double> AbsLength = {10.4 * mm};
  
  // std::vector<G4double> ScintPhoEnergy={1.95 * eV, 2.01 * eV, 2.08 * eV, 2.24 * eV, 2.51 * eV, 2.59 * eV, 2.69 * eV, 3.09 * eV, 3.40 * eV, 3.66 * eV};
  // std::vector<G4double> ScintFastArray={0.0213, 0.038, 0.0517, 0.0727, 0.0973, 0.0913, 0.078, 0.05, 0.0333, 0.02};
  // std::vector<G4double> ScintSlowArray={0.19, 0.34, 0.47, 0.65, 0.88, 0.82, 0.70, 0.45, 0.30, 0.18};


  std::vector<G4double> ScintPhoEnergy = {1.94 * eV, 1.95 * eV, 1.95 * eV, 1.96 * eV, 1.97 * eV, 1.98 * eV, 2.00 * eV, 2.01 * eV, 2.03 * eV, 2.04 * eV,
  					  2.06 * eV, 2.08 * eV, 2.10 * eV, 2.13 * eV, 2.16 * eV, 2.18 * eV, 2.21 * eV, 2.24 * eV, 2.26 * eV, 2.30 * eV,
  					  2.32 * eV, 2.36 * eV, 2.39 * eV, 2.42 * eV, 2.46 * eV, 2.51 * eV, 2.52 * eV, 2.56 * eV, 2.59 * eV, 2.63 * eV,
  					  2.65 * eV, 2.69 * eV, 2.73 * eV, 2.77 * eV, 2.81 * eV, 2.86 * eV, 2.91 * eV, 2.97 * eV, 3.04 * eV, 3.09 * eV,
  					  3.15 * eV, 3.21 * eV, 3.27 * eV, 3.33 * eV, 3.40 * eV, 3.46 * eV, 3.52 * eV, 3.57 * eV, 3.62 * eV, 3.66 * eV};
  

  std::vector<G4double> ScintFastArray = {0.019, 0.024, 0.021, 0.027, 0.029, 0.032, 0.035, 0.038, 0.042, 0.045,
  					  0.048, 0.052, 0.056, 0.060, 0.063, 0.066, 0.069, 0.073, 0.076, 0.080,
  					  0.082, 0.088, 0.090, 0.093, 0.096, 0.097, 0.096, 0.094, 0.091, 0.086,
  					  0.084, 0.078, 0.074, 0.071, 0.067, 0.064, 0.061, 0.057, 0.054, 0.050,
  					  0.046, 0.043, 0.040, 0.037, 0.033, 0.031, 0.028, 0.025, 0.023, 0.020};
  

  std::vector<G4double> ScintSlowArray = {0.171, 0.216, 0.192, 0.240, 0.264, 0.288, 0.315, 0.342, 0.375, 0.405,
  					    0.429, 0.465, 0.501, 0.537, 0.564, 0.594, 0.621, 0.654, 0.681, 0.720,
  					    0.741, 0.792, 0.807, 0.837, 0.861, 0.876, 0.867, 0.849, 0.822, 0.774,
  					    0.756, 0.702, 0.669, 0.639, 0.603, 0.573, 0.552, 0.516, 0.483, 0.450,
  					    0.417, 0.387, 0.360, 0.330, 0.300, 0.276, 0.252, 0.228, 0.204, 0.183};


  // std::vector<G4double> ScintPhoEnergy = {1.94 * eV, 2.51 * eV, 3.66 * eV};
  // std::vector<G4double> ScintFastArray = {0.019, 0.097, 0.020};
  // std::vector<G4double> ScintSlowArray = {0.171, 0.876, 0.183};

  // G4MaterialPropertyVector *MPVec;
  // MPVec->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintilSlowArray, false, true);

  G4int lenArray = 10;
  auto myMPT1 = new G4MaterialPropertiesTable();
  
  // Values can be added to the material property table individually.
  // With this method, spline interpolation cannot be set. Arguments
  // createNewKey and spline both take their default values of false.
  // Need to specify the number of entries (1) in the arrays, as an argument
  // to AddProperty.
  // G4int numEntries = 1;
  // myMPT1->AddProperty("RINDEX", &photonEnergy[0], &refractiveIndex1[0], numEntries);

  // for (size_t i = 1; i < photonEnergy.size(); ++i) {
  //   myMPT1->AddEntry("RINDEX", photonEnergy[i], refractiveIndex1[i]);
  // }

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

  // Adding a property using a C-style array.
  // Spline interpolation isn't used for scintillation.
  // Arguments spline and createNewKey both take default value false.

  //G4cout << "\n\ndatatype is: " <<typeid( myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, true)).name() << G4endl;

  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, false);
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", ScintPhoEnergy, ScintSlowArray, false, false);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 8200. / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 60. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 300. * ns); 
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.1);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.9);
  myMPT1->AddConstProperty("SCINTILLATIONRISETIME1", 30 * ps);

  G4MaterialPropertyVector *MPTVec = myMPT1->GetProperty("SCINTILLATIONCOMPONENT1");

  G4cout << "\n\nMaterial Property vector is: " << *MPTVec << G4endl;

    // Set the Birks Constant for scintillator
  //DetectorMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  //DetectorMat->GetIonisation()->SetBirksConstant(0 * mm / MeV);

  //myMPT1->DumpTable();

  
  DetectorMat->SetMaterialPropertiesTable(myMPT1);
    

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
  G4double rad_dz = 5.0 * cm;
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
 












































