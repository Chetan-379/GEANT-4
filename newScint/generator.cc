#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  G4ThreeVector pos(0.*cm, 0.*cm, -5.*cm);

  //G4ThreeVector pos(0.*cm, 0.*cm, -0.5*cm);  
  
 
   G4ThreeVector mom(0, 0, 1.);
  //G4ThreeVector mom(0, 0.5, 1.);
  

  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);  
  //fParticleGun->SetParticleMomentum(511*keV);
  fParticleGun->SetParticleEnergy(511*keV);
  //fParticleGun->SetParticleEnergy(2.5*eV);
  
    
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName="gamma";
  //G4String particleName="opticalphoton";
  G4ParticleDefinition *particle = particleTable->FindParticle(particleName);
  
  fParticleGun->SetParticleDefinition(particle);

  G4ThreeVector gammaMom = fParticleGun->GetParticleMomentumDirection();  
  // G4ThreeVector pol1 = gammaMom.orthogonal().unit();  
  // pol1 = pol1.rotate(2*CLHEP::pi * G4UniformRand(), gammaMom);     //randomizing the photon polarisation direction in the plane perp to gamma mom  
  // fParticleGun->SetParticlePolarization(pol1);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);      
}
