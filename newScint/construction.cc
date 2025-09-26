#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
 G4double energy[2] = {1.2398419*eV/0.9, 1.2398419*eV/0.2};
  G4double rindexWorld[2] = {1.0, 1.0};

  //World(Air)
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

  G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
  mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
  worldMat->SetMaterialPropertiesTable(mptWorld);
 
//Defining the GAGG scintillator (Gd3Al2Ga3O12:Ce (Cerium is just activator which affects SY and decay timing. I am not adding it in the Material explicitly))
  
  G4Material *gagg = new G4Material("GAGG", 6.63*g/cm3, 4);
  G4Element *Gd = nist->FindOrBuildElement("Gd");
  G4Element *Al = nist->FindOrBuildElement("Al");
  G4Element *Ga = nist->FindOrBuildElement("Ga");
  G4Element *O = nist->FindOrBuildElement("O");
  
  gagg -> AddElement(Gd, 3);
  gagg -> AddElement(Al, 2);
  gagg -> AddElement(Ga, 3);
  gagg -> AddElement(O, 12);

  G4Material *lyso = new G4Material("LYSO", 7.1*g/cm3, 4);  
  G4Element *Lu = nist->FindOrBuildElement("Lu");
  G4Element *Y = nist->FindOrBuildElement("Y");
  G4Element *Si = nist->FindOrBuildElement("Si");
  // G4Element *O = nist->FindOrBuildElement("O");
  
  lyso -> AddElement(Lu, 9);
  lyso -> AddElement(Y, 1);
  lyso -> AddElement(Si, 5);
  lyso -> AddElement(O, 25);

  //G4cout << "\n\nName of the Detector Material is: " << DetectorMat->GetBaseMaterial()->GetName() << G4endl;
  
 

  std::vector<G4double> G_RIPhoEnergy = {1.56 *eV, 2.19 *eV, 2.38 *eV, 3.47 *eV, 4.24 *eV, 4.65 *eV};
  std::vector<G4double> G_RefractiveIdxDet = {1.872, 1.873, 1.870, 1.887, 1.914, 1.940};

  std::vector<G4double> L_RIPhoEnergy = {0.996 *eV, 1.177 *eV, 1.425 *eV, 1.673 *eV, 1.968 *eV, 2.304 *eV, 2.643 *eV, 2.841 *eV, 2.949 *eV,
				       3.044 *eV, 3.130 *eV, 3.215 *eV, 3.328 *eV, 3.432 *eV, 3.680 *eV, 4.036 *eV, 4.383 *eV, 4.789 *eV};
  
  std::vector<G4double> L_RefractiveIdxDet = {1.784, 1.784, 1.789, 1.797, 1.805, 1.808, 1.816, 1.824, 1.833,
					    1.835, 1.833, 1.838, 1.835, 1.835, 1.841, 1.852, 1.868, 1.882};
  // std::vector<G4double> AbsPhoEnergy = {511 * keV};
  // std::vector<G4double> AbsLength = {10.4 * mm};

  std::vector<G4double> G_ScintPhoEnergy = {1.771 *eV, 1.791 *eV, 1.816 *eV, 1.834 *eV, 1.857 *eV, 1.878 *eV, 1.901 *eV, 1.916 *eV, 1.933 *eV, 1.951 *eV,
					  1.970 *eV, 1.988 *eV, 2.010 *eV, 2.025 *eV, 2.040 *eV, 2.056 *eV, 2.074 *eV, 2.091 *eV, 2.106 *eV, 2.122 *eV,
					  2.133 *eV, 2.143 *eV, 2.155 *eV, 2.164 *eV, 2.174 *eV, 2.181 *eV, 2.190 *eV, 2.196 *eV, 2.207 *eV, 2.217 *eV,
					  2.230 *eV, 2.243 *eV, 2.254 *eV, 2.264 *eV, 2.274 *eV, 2.286 *eV, 2.296 *eV, 2.310 *eV, 2.319 *eV, 2.329 *eV,
					  2.342 *eV, 2.354 *eV, 2.372 *eV, 2.389 *eV, 2.405 *eV, 2.421 *eV, 2.431 *eV, 2.444 *eV, 2.450 *eV, 2.453 *eV,
					  2.460 *eV, 2.467 *eV, 2.471 *eV, 2.475 *eV, 2.479 *eV, 2.485 *eV, 2.491 *eV, 2.496 *eV, 2.499 *eV, 2.506 *eV,
					  2.510 *eV, 2.515 *eV, 2.519 *eV, 2.522 *eV, 2.526 *eV, 2.531 *eV, 2.536 *eV, 2.539 *eV, 2.543 *eV, 2.548 *eV,
					  2.557 *eV, 2.565 *eV, 2.571 *eV, 2.580 *eV, 2.589 *eV, 2.598 *eV, 2.613 *eV, 2.631 *eV, 2.645 *eV, 2.665 *eV};

  std::vector<G4double> G_ScintFastArray = {0.032, 0.036, 0.039, 0.045, 0.051, 0.057, 0.063, 0.069, 0.073, 0.078,
					    0.086, 0.093, 0.104, 0.112, 0.122, 0.128, 0.140, 0.157, 0.174, 0.191,
					    0.201, 0.214, 0.225, 0.236, 0.244, 0.253, 0.261, 0.268, 0.278, 0.289,
					    0.301, 0.313, 0.324, 0.334, 0.343, 0.353, 0.360, 0.366, 0.374, 0.381,
					    0.385, 0.389, 0.392, 0.392, 0.391, 0.388, 0.381, 0.373, 0.364, 0.355,
					    0.345, 0.335, 0.324, 0.314, 0.306, 0.292, 0.279, 0.266, 0.254, 0.236,
					    0.218, 0.204, 0.193, 0.183, 0.170, 0.157, 0.145, 0.133, 0.121, 0.110,
					    0.099, 0.087, 0.078, 0.065, 0.051, 0.039, 0.028, 0.017, 0.008, 0.002};
  

  std::vector<G4double> G_ScintSlowArray = {0.050, 0.055, 0.061, 0.069, 0.079, 0.088, 0.097, 0.107, 0.114, 0.122,
					    0.134, 0.144, 0.161, 0.174, 0.189, 0.199, 0.217, 0.244, 0.270, 0.296,
					    0.312, 0.332, 0.348, 0.365, 0.379, 0.392, 0.405, 0.416, 0.432, 0.448,
					    0.466, 0.486, 0.502, 0.517, 0.533, 0.548, 0.558, 0.568, 0.579, 0.591,
					    0.598, 0.603, 0.608, 0.608, 0.606, 0.602, 0.591, 0.578, 0.565, 0.551,
					    0.535, 0.519, 0.503, 0.486, 0.474, 0.453, 0.433, 0.413, 0.394, 0.365,
					    0.339, 0.317, 0.300, 0.284, 0.264, 0.244, 0.224, 0.207, 0.188, 0.171,
					    0.153, 0.136, 0.121, 0.100, 0.079, 0.061, 0.044, 0.027, 0.013, 0.003};
  
  std::vector<G4double> L_ScintPhoEnergy = {2.189 *eV, 2.199 *eV, 2.212 *eV, 2.222 *eV, 2.234 *eV, 2.254 *eV, 2.268 *eV, 2.285 *eV, 2.293 *eV, 2.316 *eV,
					    2.330 *eV, 2.350 *eV, 2.366 *eV, 2.373 *eV, 2.380 *eV, 2.396 *eV, 2.408 *eV, 2.411 *eV, 2.428 *eV, 2.444 *eV,
					    2.450 *eV, 2.458 *eV, 2.471 *eV, 2.493 *eV, 2.508 *eV, 2.524 *eV, 2.534 *eV, 2.538 *eV, 2.546 *eV, 2.573 *eV,
					    2.588 *eV, 2.593 *eV, 2.603 *eV, 2.619 *eV, 2.634 *eV, 2.650 *eV, 2.657 *eV, 2.668 *eV, 2.687 *eV, 2.714 *eV,
					    2.753 *eV, 2.784 *eV, 2.815 *eV, 2.824 *eV, 2.839 *eV, 2.880 *eV, 2.914 *eV, 2.924 *eV, 2.940 *eV, 2.954 *eV,
					    2.962 *eV, 2.970 *eV, 2.984 *eV, 2.989 *eV, 3.001 *eV, 3.018 *eV, 3.028 *eV, 3.037 *eV, 3.090 *eV, 3.111 *eV,
					    3.144 *eV, 3.165 *eV, 3.195 *eV, 3.229 *eV, 3.247 *eV, 3.263 *eV, 3.285 *eV, 3.314 *eV, 3.354 *eV, 3.385 *eV,
					    3.400 *eV, 3.418 *eV, 3.435 *eV, 3.513 *eV, 3.559 *eV, 3.635 *eV};

  std::vector<G4double> L_ScintFastArray = {0.016, 0.004, 0.023, 0.002, 0.019, 0.026, 0.021, 0.019, 0.035, 0.033,
					    0.035, 0.047, 0.055, 0.039, 0.052, 0.055, 0.070, 0.064, 0.082, 0.084,
					    0.103, 0.100, 0.112, 0.127, 0.156, 0.167, 0.187, 0.204, 0.191, 0.244,
					    0.267, 0.284, 0.289, 0.323, 0.356, 0.365, 0.394, 0.411, 0.444, 0.504,
					    0.602, 0.698, 0.760, 0.791, 0.819, 0.887, 0.931, 0.941, 0.966, 0.960,
					    0.990, 0.980, 0.990, 0.996, 0.980, 0.969, 0.934, 0.929, 0.803, 0.754,
					    0.658, 0.594, 0.501, 0.372, 0.309, 0.250, 0.185, 0.106, 0.050, 0.036,
					    0.008, 0.017, 0.006, 0.016, 0.006, 0.006};

  


  auto G_MPT = new G4MaterialPropertiesTable();
  auto L_MPT = new G4MaterialPropertiesTable();
  
  G_MPT->AddProperty("RINDEX", G_RIPhoEnergy, G_RefractiveIdxDet);
  L_MPT->AddProperty("RINDEX", L_RIPhoEnergy, L_RefractiveIdxDet);
  //L_MPT->AddProperty("RINDEX", G_RIPhoEnergy, G_RefractiveIdxDet);

  // Check that group velocity is calculated from RINDEX
  if (G_MPT->GetProperty("RINDEX")->GetVectorLength()
      != G_MPT->GetProperty("GROUPVEL")->GetVectorLength())
  {
    G4ExceptionDescription ed;
    ed << "Error calculating group velocities. Incorrect number of entries "
          "in group velocity material property vector.";
    G4Exception("OpNovice::OpNoviceDetectorConstruction", "OpNovice001", FatalException, ed);
  }

  if (L_MPT->GetProperty("RINDEX")->GetVectorLength()
      != L_MPT->GetProperty("GROUPVEL")->GetVectorLength())
  {
    G4ExceptionDescription ed;
    ed << "Error calculating group velocities. Incorrect number of entries "
          "in group velocity material property vector.";
    G4Exception("OpNovice::OpNoviceDetectorConstruction", "OpNovice001", FatalException, ed);
  }

  G_MPT->AddProperty("SCINTILLATIONCOMPONENT1", G_ScintPhoEnergy, G_ScintFastArray, false, true);
  G_MPT->AddProperty("SCINTILLATIONCOMPONENT2", G_ScintPhoEnergy, G_ScintSlowArray, false, true);
  G_MPT->AddConstProperty("SCINTILLATIONYIELD", 28244 / MeV);
  G_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  G_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 50.1 * ns);
  G_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 321.5 * ns);
  G_MPT->AddConstProperty("SCINTILLATIONYIELD1", 0.392);
  G_MPT->AddConstProperty("SCINTILLATIONYIELD2", 0.608);
  G_MPT->AddConstProperty("SCINTILLATIONRISETIME1", 8 * ns);
  
  L_MPT->AddProperty("SCINTILLATIONCOMPONENT1", L_ScintPhoEnergy, L_ScintFastArray, false, true);
  L_MPT->AddConstProperty("SCINTILLATIONYIELD", 30000 / MeV);
  L_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  L_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40 * ns);
  L_MPT->AddConstProperty("SCINTILLATIONYIELD1", 1);
  L_MPT->AddConstProperty("SCINTILLATIONRISETIME1", 70 * ps);



  // Set the Birks Constant for scintillator
  gagg->GetIonisation()->SetBirksConstant(0. * mm / MeV);
  lyso->GetIonisation()->SetBirksConstant(0. * mm / MeV);
  //==============================================
  gagg->SetMaterialPropertiesTable(G_MPT);
  lyso->SetMaterialPropertiesTable(L_MPT);
  //==============================================
     

  //----------Volumes-------------
  //
  // The world
  G4Box *solidWorld = new G4Box("solidWorld", 0.25*m, 0.25*m, 0.25*m);
  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");
  G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicWorld, "physWorld", 0, false, 0, true);

  
  //Trapezoidal enclosure tiles                                                          
  G4double dim_X = 3 * mm, dim_Y = 3 * mm;  
  G4double dim_Z = 10. * mm;
   
  G4Box *solidDetector = new G4Box("solidDetector", dim_X/2, dim_Y/2, dim_Z/2); 	      

  logic_gagg = new G4LogicalVolume(solidDetector, gagg, "gagg_LV"); 
  G4VPhysicalVolume *phys_gagg = new G4PVPlacement(0, G4ThreeVector(0,0,-5*mm), logic_gagg, "GAGG_PV", logicWorld, false, 1, true);
  G4VisAttributes* visAttr_gagg = new G4VisAttributes(G4Colour(1.3, 1.2, 0.6,0.5));

  visAttr_gagg->SetForceSolid(true);     
  visAttr_gagg->SetVisibility(true);
  logic_gagg->SetVisAttributes(visAttr_gagg);


  logic_lyso = new G4LogicalVolume(solidDetector, lyso, "lyso_LV");          
  G4VPhysicalVolume *phys_lyso = new G4PVPlacement(0, G4ThreeVector(0,0,5*mm), logic_lyso, "LYSO_PV", logicWorld, false, 2, true);

  G4VisAttributes* visAttr_lyso = new G4VisAttributes(G4Colour(0.53, 0.51, 0.93,0.5));
  visAttr_lyso->SetForceSolid(true);
  visAttr_lyso->SetVisibility(true);  
  logic_lyso->SetVisAttributes(visAttr_lyso);

  
  //defining the surfaces
  G4OpticalSurface* OpSurface_ScintAir = new G4OpticalSurface("Surface_scint_air");                
  OpSurface_ScintAir->SetType(dielectric_metal);
  OpSurface_ScintAir->SetModel(unified);
  OpSurface_ScintAir->SetFinish(polished);
  
  //OpSurface->SetSigmaAlpha(0.1);
  
  std::vector<G4double> pp = {1.5 * eV, 3.7 * eV};
  //std::vector<G4double> Pho_Energy_LASurface = {2.0 * eV, 3.7 * eV};

  std::vector<G4double> rindex = {1.35, 1.40};
  
  std::vector<G4double> reflectivity = {0.99, 0.99};
  //std::vector<G4double> reflectivity = {1, 1};
  //std::vector<G4double> efficiency = {0.5, 0.5};  

  G4MaterialPropertiesTable* SMPT_ScintAir = new G4MaterialPropertiesTable();
  
  SMPT_ScintAir->AddProperty("RINDEX", pp, rindex);
  // SMPT->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe);
  // SMPT->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike);
  // SMPT->AddProperty("BACKSCATTERCONSTANT", pp, backscatter);
  // SMPT->AddProperty("REFLECTIVITY", pp, reflectivity);
  SMPT_ScintAir->AddProperty("REFLECTIVITY", pp, reflectivity);
  //SMPT_ScintAir->AddProperty("EFFICIENCY", pp, efficiency);
  
  OpSurface_ScintAir->SetMaterialPropertiesTable(SMPT_ScintAir);

 G4OpticalSurface* OpSurface_ScintScint = new G4OpticalSurface("Surface_scint_scint");
 OpSurface_ScintScint->SetType(dielectric_dielectric);
 OpSurface_ScintScint->SetModel(unified);
 OpSurface_ScintScint->SetFinish(polished);
   
 G4MaterialPropertiesTable* SMPT_ScintScint = new G4MaterialPropertiesTable();
 SMPT_ScintScint->AddProperty("RINDEX", pp, rindex);
 SMPT_ScintScint->AddProperty("REFLECTIVITY", pp, reflectivity);
 //SMPT_ScintScint->AddProperty("EFFICIENCY", pp, efficiency);

 OpSurface_ScintScint->SetMaterialPropertiesTable(SMPT_ScintScint); 
  
 G4LogicalBorderSurface *gagg_AirBoundary = new G4LogicalBorderSurface("GAGG_Air_Boundary", phys_gagg, physWorld, OpSurface_ScintAir);
 G4LogicalBorderSurface *lyso_AirBoundary = new G4LogicalBorderSurface("LYSO_Air_Boundary", phys_lyso, physWorld, OpSurface_ScintAir);
 
 //G4LogicalBorderSurface *gagg_lysoBoundary = new G4LogicalBorderSurface("GAGG_LYSO_Boundary", phys_gagg, phys_lyso, OpSurface_ScintScint);
 
  //fScoringVolume = logicDetector;
  
  return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
  // auto InStripSD = new CellSD("InStripSD", "InStripHitsCollection", 312);  //sensitive detector for the innermost layer
  // G4SDManager::GetSDMpointer()->AddNewDetector(InStripSD);
  // SetSensitiveDetector("InStripLV", InStripSD);
  // InStripSD->SetVerboseLevel(2);
  

  // auto StripSD = new CellSD("StripSD", "StripHitsCollection", 96);
  // G4SDManager::GetSDMpointer()->AddNewDetector(StripSD);
  // SetSensitiveDetector("StripLV", StripSD);
  // StripSD->SetVerboseLevel(2);

  // auto OuterMostStripSD = new CellSD("OuterMostStripSD", "OuterMostStripHitsCollection", 96);
  // G4SDManager::GetSDMpointer()->AddNewDetector(OuterMostStripSD);
  // SetSensitiveDetector("OuterMostStripLV", OuterMostStripSD);
  // OuterMostStripSD->SetVerboseLevel(2);

  // // //setting the logical volume to be unpolarised
  // G4PolarizationManager* polMgr = G4PolarizationManager::GetInstance();
  // G4cout << "\n\n ========= volpol" << polMgr->GetVolumePolarization (logicOuterMostStrip) << G4endl;
  // polMgr->SetVolumePolarization(logicInStrip, G4ThreeVector(0., 0., 0.));
  // polMgr->SetVolumePolarization(logicOuterMostStrip, G4ThreeVector(0., 0., 0.));


  
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
 












































