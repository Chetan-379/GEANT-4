#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
  fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Construction");
  fMessenger->DeclareProperty("TileDepth", depth, "Thickness of Scintillator tile");

  depth = 5;

  DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

void MyDetectorConstruction::DefineMaterials()
{
    //-----------------Defining and Adding the Material Properties----------
  //
  
  G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexWorld[2] = {1.0, 1.0};

  //World(Air)
  G4NistManager *nist = G4NistManager::Instance();
  worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);
   
  DetectorMat = new G4Material("LaBr3", 5.29*g/cm3, 2);
  
  La = nist->FindOrBuildElement("La");
  Br = nist->FindOrBuildElement("Br");
  
  DetectorMat -> AddElement(La, 1);
  DetectorMat -> AddElement(Br, 3);
  
  std::vector<G4double> RIPhoEnergy = {1.772 *eV, 1.908 *eV, 2.066 *eV, 2.255 *eV, 2.481 *eV, 2.757 *eV, 2.919 *eV, 3.102 *eV, 3.305 *eV, 3.545 *eV};
  std::vector<G4double> RefractiveIdxDet = {2.025, 2.039, 2.057, 2.080, 2.113, 2.162, 2.195, 2.237, 2.291, 2.371};

  std::vector<G4double> ScintPhoEnergy = {2.813 *eV, 2.863 *eV, 2.891 *eV, 2.919 *eV, 2.948 *eV, 2.982 *eV, 2.997 *eV, 3.010 *eV, 3.020 *eV, 3.028 *eV,
					  3.038 *eV, 3.044 *eV, 3.049 *eV, 3.060 *eV, 3.067 *eV, 3.077 *eV, 3.085 *eV, 3.093 *eV, 3.103 *eV, 3.113 *eV,
					  3.121 *eV, 3.131 *eV, 3.145 *eV, 3.160 *eV, 3.188 *eV, 3.236 *eV, 3.272 *eV, 3.302 *eV, 3.350 *eV, 3.400 *eV,
					  3.436 *eV, 3.483 *eV, 3.504 *eV, 3.520 *eV, 3.529 *eV, 3.531 *eV, 3.536 *eV, 3.544 *eV, 3.553 *eV, 3.559 *eV,
					  3.563 *eV, 3.574 *eV, 3.587 *eV, 3.596 *eV, 3.613 *eV, 3.635 *eV, 3.663 *eV, 3.713 *eV, 3.829 *eV, 3.906 *eV};

  std::vector<G4double> ScintFastArray = {0.002, 0.023, 0.046, 0.090, 0.131, 0.183, 0.212, 0.246, 0.282, 0.327,
					  0.376, 0.409, 0.448, 0.481, 0.508, 0.546, 0.582, 0.628, 0.665, 0.731,
					  0.763, 0.790, 0.819, 0.853, 0.878, 0.890, 0.879, 0.886, 0.898, 0.885,
					  0.871, 0.840, 0.798, 0.763, 0.731, 0.699, 0.649, 0.591, 0.531, 0.466,
					  0.412, 0.333, 0.298, 0.258, 0.198, 0.145, 0.094, 0.045, 0.018, 0.002};
  

  std::vector<G4double> ScintSlowArray = {0.000, 0.003, 0.005, 0.010, 0.015, 0.020, 0.024, 0.027, 0.031, 0.036,
					  0.042, 0.045, 0.050, 0.053, 0.056, 0.061, 0.065, 0.070, 0.074, 0.081,
					  0.085, 0.088, 0.091, 0.095, 0.098, 0.099, 0.098, 0.098, 0.100, 0.098,
					  0.097, 0.093, 0.089, 0.085, 0.081, 0.078, 0.072, 0.066, 0.059, 0.052,
					  0.046, 0.037, 0.033, 0.029, 0.022, 0.016, 0.010, 0.005, 0.002, 0.000};


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
  
  //myMPT1->AddProperty("ABSLENGTH", AbsPhoEnergy, AbsLength, false, true);

  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", ScintPhoEnergy, ScintFastArray, false, true);
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", ScintPhoEnergy, ScintSlowArray, false, true);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 61000 / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 30. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 300. * ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.9);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.1);
  myMPT1->AddConstProperty("SCINTILLATIONRISETIME1", 5.7 * ns);

  // Set the Birks Constant for scintillator
  DetectorMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
  //==============================================
  DetectorMat->SetMaterialPropertiesTable(myMPT1);
  //==============================================

}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{

   

  //----------Volumes-------------
  //
  // The world
  solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);
  logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");
  physWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld, "physWorld", 0, false, 0, true);
    
  //Trapezoidal enclosure tiles                                                          
  // G4double rad_dxa = 20 * cm, rad_dxb = 20 * cm;  
  //G4double rad_dz = 10 * cm;

  G4double rad_dxa = 2.5 * cm, rad_dxb = 2.5 * cm;  
  G4double rad_dz = depth *cm;
  G4double rad_dya = (rad_dxa+2*rad_dz), rad_dyb = (rad_dxa + 2 * rad_dz);

  G4RotationMatrix new_rotmX, new_rotmY, new_rotmZ;
  G4int Pos[2] = {-1,1};
  new_rotmX.rotateX(90 * deg);
  new_rotmY.rotateY(90 * deg);  
  new_rotmZ.rotateZ(90 * deg);
  
  G4ThreeVector position;
  
  solidDetector = new G4Trd("solidDetector",  // its name
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
      
      physDetector = new G4PVPlacement(transform, logicDetector, "physDetector", logicWorld, false, (i+1)*(j+1), true);

    }
  }

  // The enclosure
  solidEncVol = new G4Box("EncVol", rad_dxa/2, rad_dxa/2, rad_dxa/2);
  logicEncVol = new G4LogicalVolume(solidEncVol, worldMat, "logicalEncVol");
  physEncVol = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicEncVol, "physEncVol", logicWorld, false, 0, true);

  
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
 












































