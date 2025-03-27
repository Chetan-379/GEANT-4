#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  //-----------------Defining the and Adding the Material Properties----------
  //
  
  G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexAerogel[2] = {1.1, 1.1};
  G4double rindexWorld[2] = {1.0, 1.0};

  std::vector<G4double> photonEnergy = {
   2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV, 2.256 * eV, 2.298 * eV,
   2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV, 2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV,
   2.757 * eV, 2.820 * eV, 2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
   3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV, 4.002 * eV, 4.136 * eV};


  //World(Air)
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);

  //Detector Tile
  //G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_PbWO4");
  G4Material *DetectorMat = nist->FindOrBuildMaterial("G4_Si");
  

    std::vector<G4double> refractiveIndex1 = {
    1.3435, 1.344, 1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,  1.3475, 1.348,  1.3485,
    1.3492, 1.35,  1.3505, 1.351,  1.3518, 1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,
    1.3555, 1.356, 1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608};
  std::vector<G4double> absorption = {
    3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m, 15.152 * m, 17.241 * m,
    18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m, 45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m,
    55.556 * m, 52.632 * m, 52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
    30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m, 17.500 * m, 14.500 * m};

  // Material properties can be added as arrays. However, in this case it is
  // up to the user to make sure both arrays have the same number of elements.
  G4double scintilFastArray[]{1.0, 1.0};
  G4double energyArray[]{2.034 * eV, 4.136 * eV};
  G4int lenArray = 2;

  std::vector<G4double> scintilSlow = {
    0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00, 3.00, 2.00,
    1.00, 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 5.00, 4.00};
  
  auto myMPT1 = new G4MaterialPropertiesTable();
  
  // Values can be added to the material property table individually.
  // With this method, spline interpolation cannot be set. Arguments
  // createNewKey and spline both take their default values of false.
  // Need to specify the number of entries (1) in the arrays, as an argument
  // to AddProperty.
  G4int numEntries = 1;
  myMPT1->AddProperty("RINDEX", &photonEnergy[0], &refractiveIndex1[0], numEntries);

  for (size_t i = 1; i < photonEnergy.size(); ++i) {
    myMPT1->AddEntry("RINDEX", photonEnergy[i], refractiveIndex1[i]);
  }

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
  myMPT1->AddProperty("ABSLENGTH", photonEnergy, absorption, false, true);

  // Adding a property using a C-style array.
  // Spline interpolation isn't used for scintillation.
  // Arguments spline and createNewKey both take default value false.
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", energyArray, scintilFastArray, lenArray);
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintilSlow, false, true);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 50. / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.8);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.2);

  // std::vector<G4double> energy_water = {
  //   1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV, 1.67567 * eV,
  //   1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV, 1.7971 * eV,  1.82352 * eV,
  //   1.85074 * eV, 1.87878 * eV, 1.90769 * eV, 1.93749 * eV, 1.96825 * eV, 1.99999 * eV,
  //   2.03278 * eV, 2.06666 * eV, 2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV,
  //   2.25454 * eV, 2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
  //   2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV, 2.81817 * eV,
  //   2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV, 3.17948 * eV, 3.26315 * eV,
  //   3.35134 * eV, 3.44444 * eV, 3.54285 * eV, 3.64705 * eV, 3.75757 * eV, 3.87499 * eV,
  //   3.99999 * eV, 4.13332 * eV, 4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV,
  //   4.95999 * eV, 5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV};

  // // Rayleigh scattering length is calculated by G4OpRayleigh

  // // Mie: assume 100 times larger than the rayleigh scattering
  // std::vector<G4double> mie_water = {
  //   167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m, 128587.4 * m,
  //   121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m, 97200.35 * m, 91686.86 * m,
  //   86411.33 * m, 81366.79 * m, 76546.42 * m, 71943.46 * m, 67551.29 * m, 63363.36 * m,
  //   59373.25 * m, 55574.61 * m, 51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m,
  //   39239.39 * m, 36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
  //   24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m, 16072.45 * m,
  //   14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m, 9920.416 * m, 8941.407 * m,
  //   8036.711 * m, 7202.470 * m, 6434.927 * m, 5730.429 * m, 5085.425 * m, 4496.467 * m,
  //   3960.210 * m, 3473.413 * m, 3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m,
  //   1675.064 * m, 1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m};

  // // Mie: gforward, gbackward, forward backward ratio
  // G4double mie_water_const[3] = {0.99, 0.99, 0.8};
  
  // myMPT1->AddProperty("MIEHG", energy_water, mie_water, false, true);
  // myMPT1->AddConstProperty("MIEHG_FORWARD", mie_water_const[0]);
  // myMPT1->AddConstProperty("MIEHG_BACKWARD", mie_water_const[1]);
  // myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO", mie_water_const[2]);

  G4cout << "Water G4MaterialPropertiesTable:" << G4endl;
  myMPT1->DumpTable();
  
  DetectorMat->SetMaterialPropertiesTable(myMPT1);
  
  // Set the Birks Constant for the Water scintillator
  DetectorMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  
  //Na22 Radioactive source
  G4Isotope *Na22 = new G4Isotope("Na22", 11, 22, 21.9944 * g / mole);
  G4Element *elNa22 = new G4Element("Sodium-22", "Na22", 1);
  elNa22-> AddIsotope(Na22, 100.0 * perCent);

  G4Material *matNa22 = new G4Material("Na22Source", 0.97 * g / cm3, 1);
  matNa22->AddElement(elNa22, 100.0 * perCent);

  

  //----------Volumes-------------
  //
  // The world
  G4Box *solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");
  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld, "physWorld", 0, false, 0, true);

  // The world
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

      // //----------------Defining Surface---------------
  //     // Detectors
  //     auto opDetSurface = new G4OpticalSurface("DetSurface");
  //     opDetSurface->SetType(dielectric_LUTDAVIS);
  //     opDetSurface->SetFinish(Rough_LUT);
  //     opDetSurface->SetModel(DAVIS);

      
  //     auto DetSurface =
  //     	 new G4LogicalBorderSurface("DetSurface", physDetector, physWorld, opDetSurface);
      
  //     auto opticalSurface = dynamic_cast<G4OpticalSurface*>(
  //     DetSurface->GetSurface(physDetector, physWorld)->GetSurfaceProperty());
  //     //if (opticalSurface) opticalSurface->DumpInfo();  
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
 












































