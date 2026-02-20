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
   
//Defining the GAGG scintillator (Gd3Al2Ga3O12:Ce (Cerium is just activator which affects SY and decay timing. I am not adding it in the Material explicitly))
  
  gagg = new G4Material("GAGG", 6.63*g/cm3, 4);
  Gd = nist->FindOrBuildElement("Gd");
  Al = nist->FindOrBuildElement("Al");
  Ga = nist->FindOrBuildElement("Ga");
  O = nist->FindOrBuildElement("O");
  
  gagg -> AddElement(Gd, 3);
  gagg -> AddElement(Al, 2);
  gagg -> AddElement(Ga, 3);
  gagg -> AddElement(O, 12);

  lyso = new G4Material("LYSO", 7.1*g/cm3, 4);  
  Lu = nist->FindOrBuildElement("Lu");
  Y = nist->FindOrBuildElement("Y");
  Si = nist->FindOrBuildElement("Si");
  
  lyso -> AddElement(Lu, 9);
  lyso -> AddElement(Y, 1);
  lyso -> AddElement(Si, 5);
  lyso -> AddElement(O, 25);

  SiPM_Mat = nist->FindOrBuildMaterial("G4_Si");

  //G4cout << "\n\nName of the Detector Material is: " << DetectorMat->GetBaseMaterial()->GetName() << G4endl;
  
 

  std::vector<G4double> G_RIPhoEnergy = {1.56 *eV, 2.19 *eV, 2.38 *eV, 3.47 *eV, 4.24 *eV, 4.65 *eV};
  std::vector<G4double> G_RefractiveIdxDet = {1.872, 1.873, 1.870, 1.887, 1.914, 1.940};

  std::vector<G4double> L_RIPhoEnergy = {0.996 *eV, 1.177 *eV, 1.425 *eV, 1.673 *eV, 1.968 *eV, 2.304 *eV, 2.643 *eV, 2.841 *eV, 2.949 *eV,
				       3.044 *eV, 3.130 *eV, 3.215 *eV, 3.328 *eV, 3.432 *eV, 3.680 *eV, 4.036 *eV, 4.383 *eV, 4.789 *eV};
  
  std::vector<G4double> L_RefractiveIdxDet = {1.784, 1.784, 1.789, 1.797, 1.805, 1.808, 1.816, 1.824, 1.833,
					    1.835, 1.833, 1.838, 1.835, 1.835, 1.841, 1.852, 1.868, 1.882};
  // std::vector<G4double> AbsPhoEnergy = {511 * keV};
  // std::vector<G4double> AbsLength = {10.4 * mm};

 

  //defining the emission spectra of GAGG from csv file
  std::vector<G4double> G_ScintPhoEnergy, G_ScintFastArray, G_ScintSlowArray;
  //std::ifstream GAGG_file("../test_GAGG_emission.csv");
  std::ifstream GAGG_file("../GAGG_emission.csv"); 
  std::string line;
  
  std::getline(GAGG_file, line);  //skipping the headers in csv file

  while (std::getline(GAGG_file, line)) {
    std::stringstream ss(line);
    G4double energy, FastScint, SlowScint;
    char comma;
    if (ss >> energy >> comma >> FastScint >> comma >> SlowScint) {
      G_ScintPhoEnergy.push_back(energy*eV);
      G_ScintFastArray.push_back(FastScint);
      G_ScintSlowArray.push_back(SlowScint);            
    }
  }  
  GAGG_file.close();
  //----------
  
  //defining the emission spectra of LYSO from csv file
  std::vector<G4double> L_ScintPhoEnergy, L_ScintFastArray;
  //std::ifstream LYSO_file("../test_LYSO_emission.csv");
  std::ifstream LYSO_file("../LYSO_emission.csv");
  
  std::getline(LYSO_file, line);  //skipping the headers in csv file
  while (std::getline(LYSO_file, line)) {
    std::stringstream ss(line);
    G4double energy, FastScint;
    char comma;
    if (ss >> energy >> comma >> FastScint) {
      L_ScintPhoEnergy.push_back(energy*eV);
      L_ScintFastArray.push_back(FastScint);      
    }
  }  
  LYSO_file.close();
  //------------
  
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
  // G_MPT->AddConstProperty("SCINTILLATIONYIELD", 50 / MeV);  
  G_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  G_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 50.1 * ns);
  G_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 321.5 * ns);
  G_MPT->AddConstProperty("SCINTILLATIONYIELD1", 0.392);
  G_MPT->AddConstProperty("SCINTILLATIONYIELD2", 0.608);
  G_MPT->AddConstProperty("SCINTILLATIONRISETIME1", 8 * ns);
  
  L_MPT->AddProperty("SCINTILLATIONCOMPONENT1", L_ScintPhoEnergy, L_ScintFastArray, false, true);
  L_MPT->AddConstProperty("SCINTILLATIONYIELD", 30000 / MeV);
  // L_MPT->AddConstProperty("SCINTILLATIONYIELD", 50 / MeV);
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

  //Defining the optical surfaces for Scintillators
  OpSurface_ScintAir = new G4OpticalSurface("Surface_scint_air");                
  OpSurface_ScintAir->SetType(dielectric_metal);
  OpSurface_ScintAir->SetModel(unified);
  OpSurface_ScintAir->SetFinish(polished);
  
  std::vector<G4double> pp = {1.5 * eV, 3.7 * eV};  
  std::vector<G4double> scint_rindex = {1.35, 1.40};  
  std::vector<G4double> scint_reflectivity = {0.99, 0.99};
  std::vector<G4double> scint_EFF = {0., 0.};
  
  auto *SMPT_ScintAir = new G4MaterialPropertiesTable();  
  SMPT_ScintAir->AddProperty("RINDEX", pp, scint_rindex);
  SMPT_ScintAir->AddProperty("REFLECTIVITY", pp, scint_reflectivity);
  SMPT_ScintAir->AddProperty("EFFICIENCY", pp, scint_EFF);  
  
  OpSurface_ScintAir->SetMaterialPropertiesTable(SMPT_ScintAir);

  //Optical surface for SiPM
  
  //std::vector<G4double> SiPM_EFF = {0.3, 0.3};

  //defining the Quantum Efficiency of SiPM from csv file
  std::vector<G4double> SiPM_PhoEnergy, SiPM_EFF;
  //std::ifstream SiPM_file("../SiPM_QE.csv");
  std::ifstream SiPM_file("../Abs_frac_SiPM_QE.csv");
 
  
  std::getline(SiPM_file, line);  //skipping the headers in csv file
  while (std::getline(SiPM_file, line)) {
    std::stringstream ss(line);
    G4double energy, QEFF;
    char comma;
    if (ss >> energy >> comma >> QEFF) {
      SiPM_PhoEnergy.push_back(energy*eV);
      SiPM_EFF.push_back(QEFF);      
    }
  }  
  SiPM_file.close();
  //------------
  
  std::vector<G4double> SiPM_ReR = {1.92, 1.92};
  std::vector<G4double> SiPM_ImR = {1.69, 1.69};
  
  auto *SMPT_SiPM = new G4MaterialPropertiesTable();  
  //SMPT_SiPM->AddProperty("EFFICIENCY", pp, SiPM_EFF);
  SMPT_SiPM->AddProperty("EFFICIENCY", SiPM_PhoEnergy, SiPM_EFF);
  SMPT_SiPM->AddProperty("REALRINDEX", pp, SiPM_ReR);  
  SMPT_SiPM->AddProperty("IMAGINARYRINDEX", pp, SiPM_ImR);
  
  OpSurface_SiPM = new G4OpticalSurface("Surface_SiPM", glisur, polished, dielectric_metal);
  OpSurface_SiPM->SetMaterialPropertiesTable(SMPT_SiPM);
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  //------------World volume--------------
  G4Box* solidWorld = new G4Box("WorldS", 0.5*m, 0.5*m, 0.5*m);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World_LV");
  auto physWorld = new G4PVPlacement(0, {}, logicWorld, "World_PV", 0, false, 0);


  //------------Crystal dimensions----------
  G4double X_dim = 6*mm;
  G4double Y_dim = 6*mm;
  G4double Z_dim = 20*mm;
  G4double gapThick = 0.2*mm;

  G4int MatSize = 8;
  G4double sep = 3*cm;

  G4double layerX = X_dim + gapThick;
  G4double layerY = Y_dim + gapThick;
   
  std::vector<G4int> PosIdx = {-1, 1};

  //-------------------------creating the geomtery ------------------------
  auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/4);
  
  auto logic_gagg = new G4LogicalVolume(ScintS, gagg, "GAGG_LV");
  G4VisAttributes* visAttr_gagg = new G4VisAttributes(G4Colour(1.3, 1.2, 0.6,1.));
  visAttr_gagg->SetForceSolid(true);     
  visAttr_gagg->SetVisibility(true);
  logic_gagg->SetVisAttributes(visAttr_gagg);
  
  //-----
  auto logic_lyso = new G4LogicalVolume(ScintS, lyso, "LYSO_LV");
  //auto logic_lyso = new G4LogicalVolume(ScintS, gagg, "LYSO_LV");
  G4VisAttributes* visAttr_lyso = new G4VisAttributes(G4Colour(0.53, 0.51, 0.93,1.));
  visAttr_lyso->SetForceSolid(true);
  visAttr_lyso->SetVisibility(true);  
  logic_lyso->SetVisAttributes(visAttr_lyso);


  //GAGG---------------
  auto *phys_gagg = new G4PVPlacement(nullptr,
				      G4ThreeVector(0.,0.,2*cm),
				      logic_gagg,
				      "GAGG_PV",
				      logicWorld,
				      false,
				      111,
				      false);
  
  //LYSO---------------
  auto *phys_lyso = new G4PVPlacement(nullptr,			 
				      G4ThreeVector(0.,0.,3*cm),
				      logic_lyso,
				      "LYSO_PV",
				      logicWorld,
				      false,
				      211,
				      false);

  //--------
  //SiPMs-------------
  auto SiPM_S = new G4Box("SiPM", X_dim/2, Y_dim/2, Z_dim/20);
  auto logic_SiPM = new G4LogicalVolume(SiPM_S, SiPM_Mat, "SiPM_LV");

  auto *phys_SiPM = new G4PVPlacement(nullptr,
				      G4ThreeVector(0.,0.,3.6*cm),
				      logic_SiPM,
				      "SiPM_PV",
				      logicWorld,
				      false,
				      311,
				      false);
  

  
  //-----------Placing the reflectors------------	
  new G4LogicalBorderSurface("GAGG_Air_Boundary", phys_gagg, physWorld, OpSurface_ScintAir);
  new G4LogicalBorderSurface("LYSO_Air_Boundary", phys_lyso, physWorld, OpSurface_ScintAir);
  new G4LogicalBorderSurface("LYSO_SiPM_Boundary", phys_lyso, phys_SiPM, OpSurface_SiPM);

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
 


  // for (int i =0; i<MatSize; i++){
  //   for (int j=0; j<MatSize; j++){
  //     for (int k =0; k<2; k++){
  // 	//GAGG---------------
  // 	auto *phys_gagg = new G4PVPlacement(nullptr,
  // 					    G4ThreeVector(((-7./2)+i)*(X_dim+gapThick), ((-7./2)+j)*(X_dim+gapThick),PosIdx[k]*(sep/2+Z_dim/4)),
  // 					    logic_gagg,
  // 					    "GAGG_PV",
  // 					    logicWorld,
  // 					    false,
  // 					    PosIdx[k]*(i+1+(10*(j+1))+100),
  // 					    false);
	
  // 	//LYSO---------------
  // 	auto *phys_lyso = new G4PVPlacement(nullptr,			 
  // 					    G4ThreeVector(((-7./2)+i)*(X_dim+gapThick), ((-7./2)+j)*(X_dim+gapThick),PosIdx[k]*(sep/2+(3.*Z_dim/4))),
  // 					    logic_lyso,
  // 					    "LYSO_PV",
  // 					    logicWorld,
  // 					    false,
  // 					    PosIdx[k]*(i+1+(10*(j+1))+200),
  // 					    false);

	
  // 	//-----------Placing the reflectors------------	
  // 	auto *gagg_AirBoundary = new G4LogicalBorderSurface("GAGG_Air_Boundary", phys_gagg, physWorld, OpSurface_ScintAir);
  // 	auto *lyso_AirBoundary = new G4LogicalBorderSurface("LYSO_Air_Boundary", phys_lyso, physWorld, OpSurface_ScintAir);
	
  //     }
  //   }
  // }






























