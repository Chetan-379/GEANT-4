#include "construction.hh"
#include "G4RotationMatrix.hh"
#include <bits/stdc++.h>

MyDetectorConstruction::MyDetectorConstruction()
{}

MyDetectorConstruction::~MyDetectorConstruction()
{}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  auto nist = G4NistManager::Instance();
  auto air = nist->FindOrBuildMaterial("G4_AIR");
  auto plastic = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // World
  G4double worldSize = 1.5 * m;
  auto solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
  auto logicWorld = new G4LogicalVolume(solidWorld, air, "World");
  auto physWorld = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

  // Scintillator strip (long axis along Z)
  G4double stripLength = 500.0 * mm;
  G4double stripWidth  = 19.0 * mm;
  G4double stripThick  = 7.0 * mm;

  G4double InstripLength = 500.0 * mm;
  G4double InstripWidth  = 24.0 * mm;
  G4double InstripThick  = 6.0 * mm;
        
  auto solidStrip = new G4Box("Strip", stripWidth / 2, stripThick / 2, stripLength / 2);
  logicStrip = new G4LogicalVolume(solidStrip, plastic, "StripLV");

  auto solidInStrip = new G4Box("InStrip", InstripWidth / 2, InstripThick / 2, InstripLength / 2);
  logicInStrip = new G4LogicalVolume(solidInStrip, plastic, "InStripLV");

  G4int nstrips_block = 13;

    
  // J-PET four-layer geometry parameters
  struct Layer { G4int nStrips; G4double radius; };
  std::vector<Layer> layers = {
			       {312, 362.0 * mm},				 
			       {48, 425.0 * mm},
			       {48, 467.5 * mm},
			       {96, 575.0 * mm}
  };

  G4double phi_shift = 0;
  for (size_t i = 0; i < layers.size(); ++i) {
    auto layer = layers[i];
    G4double dPhi = 2.0 * CLHEP::pi / layer.nStrips;	
    for (G4int j = 0; j < layer.nStrips; ++j) {
      G4double phi = j * dPhi;
      G4double x = layer.radius * std::cos(phi);
      G4double y = layer.radius * std::sin(phi);
	
      G4RotationMatrix rot = G4RotationMatrix();
      rot.rotateZ(phi+phi_shift);

      G4ThreeVector uz = G4ThreeVector(std::cos(phi+phi_shift), std::sin(phi+phi_shift), 0.);
      if (i == 0){
	G4ThreeVector position = (layer.radius + 0.5 * InstripWidth) * uz;
	G4Transform3D transform = G4Transform3D(rot, position);	      
	new G4PVPlacement(transform,
			  logicInStrip,
			  "InStrip",
			  logicWorld,
			  false,
			  i * 1000 + j,
			  true);
      }

      else{
	G4ThreeVector position = (layer.radius + 0.5 * stripWidth) * uz;
	G4Transform3D transform = G4Transform3D(rot, position);	      
	new G4PVPlacement(transform,
			  logicStrip,
			  "Strip",
			  logicWorld,
			  false,
			  i * 1000 + j,
			  true);
      }	      
    }
	
    G4cout << "\n\nphi_shift: " << phi_shift << G4endl;
    if ( i==1) phi_shift =  phi_shift + (dPhi/2);
    else if (i == 2) phi_shift = phi_shift + (dPhi/3); 
  }

  //   //Reading/writing GDML
  //G4GDMLParser parser;
  //parser.Write("JPET_Geometry.gdml", physWorld);
  //parser.Read("JPET_Geometry.gdml");
  //G4VPhysicalVolume* physWorld = parser.GetWorldVolume();

  //G4LogicalVolume* logicWorld = physWorld->GetLogicalVolume();

  //logicWorld->SetVisAttributes(nullptr); //making the world volume visible

    
  return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
  auto InStripSD = new CellSD("InStripSD", "InStripHitsCollection", 312);  //sensitive detector for the innermost layer
  G4SDManager::GetSDMpointer()->AddNewDetector(InStripSD);
  SetSensitiveDetector("InStripLV", InStripSD);
  InStripSD->SetVerboseLevel(2);
  

  auto StripSD = new CellSD("StripSD", "StripHitsCollection", 192);
  G4SDManager::GetSDMpointer()->AddNewDetector(StripSD);
  SetSensitiveDetector("StripLV", StripSD);
  StripSD->SetVerboseLevel(2);
  
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
 












































