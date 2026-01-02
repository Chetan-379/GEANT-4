#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{
  DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{}

void MyDetectorConstruction::DefineMaterials()
{
  //-----------------Defining and Adding the Material Properties----------  
  G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexWorld[2] = {1.0, 1.0};

  G4NistManager *nist = G4NistManager::Instance();
  //worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);
   
  DetectorMat = new G4Material("GAGG", 6.63*g/cm3, 4);
  
  Gd = nist->FindOrBuildElement("Gd");
  Al = nist->FindOrBuildElement("Al");
  Ga = nist->FindOrBuildElement("Ga");
  O = nist->FindOrBuildElement("O");
  
  DetectorMat -> AddElement(Gd, 3);
  DetectorMat -> AddElement(Al, 2);
  DetectorMat -> AddElement(Ga, 3);
  DetectorMat -> AddElement(O, 12);   

  std::vector<G4double> RIPhoEnergy = {1.56 *eV, 2.19 *eV, 2.38 *eV, 3.47 *eV, 4.24 *eV, 4.65 *eV};
  std::vector<G4double> RefractiveIdxDet = {1.872, 1.873, 1.870, 1.887, 1.914, 1.940};

  std::vector<G4double> ScintPhoEnergy = {1.771 *eV, 1.791 *eV, 1.816 *eV, 1.834 *eV, 1.857 *eV, 1.878 *eV, 1.901 *eV, 1.916 *eV, 1.933 *eV, 1.951 *eV,
					  1.970 *eV, 1.988 *eV, 2.010 *eV, 2.025 *eV, 2.040 *eV, 2.056 *eV, 2.074 *eV, 2.091 *eV, 2.106 *eV, 2.122 *eV,
					  2.133 *eV, 2.143 *eV, 2.155 *eV, 2.164 *eV, 2.174 *eV, 2.181 *eV, 2.190 *eV, 2.196 *eV, 2.207 *eV, 2.217 *eV,
					  2.230 *eV, 2.243 *eV, 2.254 *eV, 2.264 *eV, 2.274 *eV, 2.286 *eV, 2.296 *eV, 2.310 *eV, 2.319 *eV, 2.329 *eV,
					  2.342 *eV, 2.354 *eV, 2.372 *eV, 2.389 *eV, 2.405 *eV, 2.421 *eV, 2.431 *eV, 2.444 *eV, 2.450 *eV, 2.453 *eV,
					  2.460 *eV, 2.467 *eV, 2.471 *eV, 2.475 *eV, 2.479 *eV, 2.485 *eV, 2.491 *eV, 2.496 *eV, 2.499 *eV, 2.506 *eV,
					  2.510 *eV, 2.515 *eV, 2.519 *eV, 2.522 *eV, 2.526 *eV, 2.531 *eV, 2.536 *eV, 2.539 *eV, 2.543 *eV, 2.548 *eV,
					  2.557 *eV, 2.565 *eV, 2.571 *eV, 2.580 *eV, 2.589 *eV, 2.598 *eV, 2.613 *eV, 2.631 *eV, 2.645 *eV, 2.665 *eV};

  std::vector<G4double> ScintFastArray = {0.032, 0.036, 0.039, 0.045, 0.051, 0.057, 0.063, 0.069, 0.073, 0.078,
					  0.086, 0.093, 0.104, 0.112, 0.122, 0.128, 0.140, 0.157, 0.174, 0.191,
					  0.201, 0.214, 0.225, 0.236, 0.244, 0.253, 0.261, 0.268, 0.278, 0.289,
					  0.301, 0.313, 0.324, 0.334, 0.343, 0.353, 0.360, 0.366, 0.374, 0.381,
					  0.385, 0.389, 0.392, 0.392, 0.391, 0.388, 0.381, 0.373, 0.364, 0.355,
					  0.345, 0.335, 0.324, 0.314, 0.306, 0.292, 0.279, 0.266, 0.254, 0.236,
					  0.218, 0.204, 0.193, 0.183, 0.170, 0.157, 0.145, 0.133, 0.121, 0.110,
					  0.099, 0.087, 0.078, 0.065, 0.051, 0.039, 0.028, 0.017, 0.008, 0.002};
  

  std::vector<G4double> ScintSlowArray = {0.050, 0.055, 0.061, 0.069, 0.079, 0.088, 0.097, 0.107, 0.114, 0.122,
					  0.134, 0.144, 0.161, 0.174, 0.189, 0.199, 0.217, 0.244, 0.270, 0.296,
					  0.312, 0.332, 0.348, 0.365, 0.379, 0.392, 0.405, 0.416, 0.432, 0.448,
					  0.466, 0.486, 0.502, 0.517, 0.533, 0.548, 0.558, 0.568, 0.579, 0.591,
					  0.598, 0.603, 0.608, 0.608, 0.606, 0.602, 0.591, 0.578, 0.565, 0.551,
					  0.535, 0.519, 0.503, 0.486, 0.474, 0.453, 0.433, 0.413, 0.394, 0.365,
					  0.339, 0.317, 0.300, 0.284, 0.264, 0.244, 0.224, 0.207, 0.188, 0.171,
					  0.153, 0.136, 0.121, 0.100, 0.079, 0.061, 0.044, 0.027, 0.013, 0.003};

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
  myMPT1->AddConstProperty("SCINTILLATIONYIELD", 28244 / MeV);
  //myMPT1->AddConstProperty("SCINTILLATIONYIELD", 200 / MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE", 1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 50.1 * ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 321.5 * ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1", 0.392);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD2", 0.608);
  myMPT1->AddConstProperty("SCINTILLATIONRISETIME1", 8 * ns);

  // Set the Birks Constant for scintillator
  DetectorMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);

  //==============================================
  //DetectorMat->SetMaterialPropertiesTable(myMPT1);
  //==============================================
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  //------------World volume--------------
  G4Box* solidWorld = new G4Box("WorldS", 0.5*m, 0.5*m, 0.5*m);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World_LV");
  auto physWorld = new G4PVPlacement(0, {}, logicWorld, "World_PV", 0, false, 0);


  //------------Crystal dimensions----------
  G4double X_dim = 3*mm;
  G4double Y_dim = 3*mm;
  G4double Z_dim = 20*mm;
  G4double gapThick = 0.2*mm;

  G4int MatSize = 8;
  G4double sep = 3*cm;

  G4double layerX = X_dim + gapThick;
  G4double layerY = Y_dim + gapThick;
   
  std::vector<G4int> PosIdx = {-1, 1};

 //--------------------------Working Don't edit----------------------------------------------------------------------------------
  //creating the geomtery using for loop
  auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/2);
  auto ScintLV = new G4LogicalVolume(ScintS, DetectorMat, "Scint_LV");

  for (int i =0; i<MatSize; i++){
    for (int j=0; j<MatSize; j++){
  // for (int i =0; i<1; i++){
  //   for (int j=0; j<1; j++){ 

      for (int k =0; k<2; k++){
      // auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/2);
      // auto ScintLV = new G4LogicalVolume(ScintS, DetectorMat, "Scint_LV");
      
      // new G4PVPlacement(nullptr,
      // 			G4ThreeVector(((-7/2)+i)*(X_dim+gapThick), ((-7/2)+i)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
      // 			ScintLV,
      // 			"Scint_PV",
      // 			logicWorld,
      // 			false,
      // 			PosIdx[k]*(i+1+(10*(j+1))),
      // 			true);

	//G4cout << "\n ###############gap+crys: " << (X_dim+gapThick)*(0.5) <<"----------------------" << G4endl;
	//G4cout << "\n ###############gap+crys: " << (1./2) <<"----------------------" << G4endl;
      new G4PVPlacement(nullptr,
  			G4ThreeVector(((-7./2)+i)*(X_dim+gapThick), ((-7./2)+j)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
			//G4ThreeVector((X_dim+gapThick)*(1/2), ((-1/2)+j)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
			//G4ThreeVector(-7*(X_dim+gapThick)/2+i*(X_dim+gapThick),-7*(X_dim+gapThick)/2+j*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),			
  			ScintLV,
  			"Scint_PV",
  			logicWorld,
  			false,
  			PosIdx[k]*(i+1+(10*(j+1))),
  			false);

      // new G4PVPlacement(nullptr,
      // 			G4ThreeVector(-1*((2*i+1)/2)*(X_dim+gapThick), -1*((2*j+1)/2)*(Y_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
      // 			ScintLV,
      // 			"Scint_PV",
      // 			logicWorld,
      // 			false,
      // 			PosIdx[k]*(i+1+(10*(j+1))),
      // 			true);

      }
    }
  }

  //
  //-----------Placing the reflectors------------
  //
  G4OpticalSurface* OpSurface = new G4OpticalSurface("Surface_scint_air");
  OpSurface->SetType(dielectric_metal);
  OpSurface->SetModel(unified);
  OpSurface->SetFinish(polished);
  
  std::vector<G4double> pp = {1.7*eV, 1.8*eV, 2.0*eV, 2.4*eV, 2.8*eV, 3*eV};  
  std::vector<G4double> rindex = {1.35, 1.40};    
  std::vector<G4double> reflectivity = {0.99, 0.99, 0.99, 0.99, 0.99, 0.99};   
  
  G4MaterialPropertiesTable* Surf_MPT = new G4MaterialPropertiesTable();    
  //Surf_MPT->AddProperty("RINDEX", pp, rindex);                                                                                                                                      
  Surf_MPT->AddProperty("REFLECTIVITY", pp, reflectivity);
  
  OpSurface->SetMaterialPropertiesTable(Surf_MPT);
  
  G4LogicalSkinSurface * ScintSurface = new G4LogicalSkinSurface("reflector",ScintLV,OpSurface);
  
  G4VisAttributes* VisAttr = new G4VisAttributes(G4Colour(1., 1., 0., 0.5));
  VisAttr->SetForceSolid(true);   
  //ScintLV->SetVisAttributes(VisAttr);
  
  
  return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
  // auto ScintSD = new CellSD("ScintSD", "ScintHitsCollection", 128);  //sensitive detector for the scintillator crystals
  // G4SDManager::GetSDMpointer()->AddNewDetector(ScintSD);
  // //SetSensitiveDetector("Scint_LV", ScintSD);
  // SetSensitiveDetector("World_LV", ScintSD);
  // ScintSD->SetVerboseLevel(2);

  
  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  // G4ThreeVector fieldValue;
  // fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  // fMagFieldMessenger->SetVerboseLevel(1);

  // // Register the field messenger for deleting
  // G4AutoDelete::Register(fMagFieldMessenger);
}

  // new G4PVPlacement(nullptr,
  // 		       G4ThreeVector(-(layerX/2 - X_dim/2) , (layerY/2 - Y_dim/2), 0.),
  // 		       ScintLV,
  // 		       "Scint_PV",
  // 		       logicY,
  // 		       false,
  // 		       gapPV->GetCopyNo(),
  // 		       true);

















//==================BACK UP===================================
 //--------------------------Working Don't edit----------------------------------------------------------------------------------
  // //creating the geomtery using for loop
  // auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/2);
  // auto ScintLV = new G4LogicalVolume(ScintS, DetectorMat, "Scint_LV");

  // for (int i =0; i<MatSize; i++){
  //   for (int j=0; j<MatSize; j++){
  //     for (int k =0; k<2; k++){
  //     // auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/2);
  //     // auto ScintLV = new G4LogicalVolume(ScintS, DetectorMat, "Scint_LV");
      
  //     // new G4PVPlacement(nullptr,
  //     // 			G4ThreeVector(((-7/2)+i)*(X_dim+gapThick), ((-7/2)+i)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
  //     // 			ScintLV,
  //     // 			"Scint_PV",
  //     // 			logicWorld,
  //     // 			false,
  //     // 			PosIdx[k]*(i+1+(10*(j+1))),
  //     // 			true);

  //     new G4PVPlacement(nullptr,
  // 			G4ThreeVector(((-7/2)+i)*(X_dim+gapThick), ((7/2)-j)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
  // 			ScintLV,
  // 			"Scint_PV",
  // 			logicWorld,
  // 			false,
  // 			PosIdx[k]*(i+1+(10*(j+1))),
  // 			true);

  //     // new G4PVPlacement(nullptr,
  //     // 			G4ThreeVector(-1*((2*i+1)/2)*(X_dim+gapThick), -1*((2*j+1)/2)*(Y_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/2)),
  //     // 			ScintLV,
  //     // 			"Scint_PV",
  //     // 			logicWorld,
  //     // 			false,
  //     // 			PosIdx[k]*(i+1+(10*(j+1))),
  //     // 			true);

  //     // //
  //     // //-----------Placing the reflectors------------
  //     // //
  //     // G4OpticalSurface* OpSurface = new G4OpticalSurface("Surface_scint_air");
  //     // OpSurface->SetType(dielectric_metal);
  //     // OpSurface->SetModel(unified);
  //     // OpSurface->SetFinish(polished);
      
  //     // std::vector<G4double> pp = {1.7*eV, 1.8*eV, 2.0*eV, 2.4*eV, 2.8*eV, 3*eV};  
  //     // std::vector<G4double> rindex = {1.35, 1.40};    
  //     // std::vector<G4double> reflectivity = {0.99, 0.99, 0.99, 0.99, 0.99, 0.99};   
      
  //     // G4MaterialPropertiesTable* Surf_MPT = new G4MaterialPropertiesTable();    
  //     // //Surf_MPT->AddProperty("RINDEX", pp, rindex);                                                                                                                                      
  //     // Surf_MPT->AddProperty("REFLECTIVITY", pp, reflectivity);
      
  //     // OpSurface->SetMaterialPropertiesTable(Surf_MPT);
      
  //     // G4LogicalSkinSurface * ScintSurface = new G4LogicalSkinSurface("reflector",ScintLV,OpSurface);
  //     // new G4LogicalBorderSurface("reflector1", ScintPV, gapPV, OpSurface);
  //     // new G4LogicalBorderSurface("reflector2", ScintPV, physWorld, OpSurface);
  //     }
  //   }
  // }

  //-----------------------------------------------------------------------------------------------------------------------------------------------------   
 

 












































