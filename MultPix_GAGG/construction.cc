#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  // Material
   auto nist = G4NistManager::Instance();
   G4Material* scint = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
   G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

   // Mother volume (world)
   G4Box* solidWorld = new G4Box("World", 0.5*m, 0.5*m, 0.5*m);
   G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, air, "World");
   auto physWorld = new G4PVPlacement(0, {}, logicWorld, "World", 0, false, 0);

   auto Scint_mat = new G4Material("GAGG", 6.63*g/cm3, 4);
  
   auto Gd = nist->FindOrBuildElement("Gd");
   auto Al = nist->FindOrBuildElement("Al");
   auto Ga = nist->FindOrBuildElement("Ga");
   auto O = nist->FindOrBuildElement("O");
  
   Scint_mat -> AddElement(Gd, 3);
   Scint_mat -> AddElement(Al, 2);
   Scint_mat -> AddElement(Ga, 3);
   Scint_mat -> AddElement(O, 12);   


   G4double X_dim = 3*mm;
   G4double Y_dim = 3*mm;
   G4double Z_dim = 20*mm;
   G4double gapThick = 0.2*mm;

   G4int MatSize = 8;
   G4double sep = 3*cm;

   G4double layerX = X_dim + gapThick;
   G4double layerY = Y_dim + gapThick;
   
   std::vector<G4int> PosIdx = {-1, 1};

   for (int i =0; i<2; i++){     
     // Big cube (mother for replicas)
     G4Box* solidBlock = new G4Box("Block", layerX*MatSize/2, layerY*MatSize/2, Z_dim/2);
     G4LogicalVolume* logicBlock = new G4LogicalVolume(solidBlock, Scint_mat, "Block_LV");
     new G4PVPlacement(0, G4ThreeVector(0,0,PosIdx[i]*((sep/2)+(Z_dim/2))), logicBlock, "Block_PV", logicWorld, false, 0);
   
     //replicate along X
   G4Box* solidX = new G4Box("Xslice", layerX/2, layerY*MatSize/2, Z_dim/2);
   G4LogicalVolume* logicX = new G4LogicalVolume(solidX, Scint_mat, "Xslice");
   new G4PVReplica("Xrep", logicX, logicBlock, kXAxis, MatSize, layerX);
   
   //replicate along Y
   G4Box* solidY = new G4Box("crystal", layerX/2, layerY/2, Z_dim/2);
   G4LogicalVolume* logicY = new G4LogicalVolume(solidY, Scint_mat, "Scint_LV");
   new G4PVReplica("Scint_Cryst", logicY, logicX, kYAxis, MatSize, layerY);

   
   //  // // Step 3: replicate along Z
   // G4Box* solidZ = new G4Box("Zslice", 1*cm, 1*cm, 5*cm);
   // G4LogicalVolume* logicCrystal = new G4LogicalVolume(solidZ, scint, "Crystal");
   // new G4PVReplica("Zrep", logicCrystal, logicY, kZAxis, 2, 5*cm);

  //-------------placing the crystals in the layers
   auto ScintS = new G4Box("Scint", X_dim/2, Y_dim/2, Z_dim/2);  // its size

   auto ScintLV = new G4LogicalVolume(ScintS,  // its solid
   				    scint,  // its material
   				    "Scint_LV");  // its name
   
   auto XgapPV = new G4PVPlacement(nullptr,  // no rotation
   				  G4ThreeVector(-(layerX/2 - X_dim/2) , (layerY/2 - Y_dim/2), 0.),  // its position
   				  ScintLV,  // its logical volume
   			      "Scint_PV",  // its name
   			      logicY,  // its mother  volume
   			      false,  // no boolean operation
   			      0,  // copy number
   			      true);  // checking overlaps

   
  //  auto YgapS = new G4Box("YGap", X_dim/2, gapThick/2, Z_dim/2);  // its size

  //  auto YgapLV = new G4LogicalVolume(YgapS,  // its solid
  // 				    air,  // its material
  // 				    "YGapLV");  // its name
   
  //  auto YgapPV = new G4PVPlacement(nullptr,  // no rotation
  // 				   G4ThreeVector(0 , -(Y_dim/2 - gapThick/2), 0.),  // its position
  // 				  YgapLV,  // its logical volume
  // 			      "YGapPV",  // its name
  // 			      logicY,  // its mother  volume
  // 			      false,  // no boolean operation
  // 			      0,  // copy number
  // 			      false);  // checking overlaps


   // //-------------placing the gap between the crystals
  //  auto XgapS = new G4Box("XGap", gapThick/2, Y_dim/2, Z_dim/2);  // its size

  //  auto XgapLV = new G4LogicalVolume(XgapS,  // its solid
  //  				    air,  // its material
  //  				    "XGapLV");  // its name
   
  //  auto XgapPV = new G4PVPlacement(nullptr,  // no rotation
  //  				  G4ThreeVector(-(X_dim/2 - gapThick/2) , 0., 0.),  // its position
  //  				  XgapLV,  // its logical volume
  //  			      "XGapPV",  // its name
  //  			      logicY,  // its mother  volume
  //  			      false,  // no boolean operation
  //  			      0,  // copy number
  //  			      true);  // checking overlaps

   
  //  auto YgapS = new G4Box("YGap", X_dim/2, gapThick/2, Z_dim/2);  // its size

  //  auto YgapLV = new G4LogicalVolume(YgapS,  // its solid
  // 				    air,  // its material
  // 				    "YGapLV");  // its name
   
  //  auto YgapPV = new G4PVPlacement(nullptr,  // no rotation
  // 				   G4ThreeVector(0 , -(Y_dim/2 - gapThick/2), 0.),  // its position
  // 				  YgapLV,  // its logical volume
  // 			      "YGapPV",  // its name
  // 			      logicY,  // its mother  volume
  // 			      false,  // no boolean operation
  // 			      0,  // copy number
  // 			      false);  // checking overlaps

  //  //XgapLV->SetVisAttributes(G4VisAttributes::GetInvisible());
   logicBlock->SetVisAttributes(G4VisAttributes::GetInvisible());
   logicX->SetVisAttributes(G4VisAttributes::GetInvisible());
   logicY->SetVisAttributes(G4VisAttributes::GetInvisible());

   G4VisAttributes* VisAttr = new G4VisAttributes(G4Colour(0., 1., 1.));
   VisAttr->SetForceSolid(true);   
   // ScintLV->SetVisAttributes(VisAttr);
   //YgapLV->SetVisAttributes(VisAttr);
    }

   
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
 












































