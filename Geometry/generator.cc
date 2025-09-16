#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  gParticleGun = new G4ParticleGun(1);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
  delete gParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);  
  
  //G4ThreeVector mom(1.5, 2.5, 0.);
  G4ThreeVector mom(0, 1, 0.);
  G4ThreeVector ZAxis(0,0,1);
  G4ThreeVector XAxis(1,0,0);
  
  mom = mom.rotate(0.6044*((2*G4UniformRand()) - 1.0), XAxis);   //first randomizing in z
  mom = mom.rotate(2*CLHEP::pi * G4UniformRand(), ZAxis);        //then randomizing in xy plane   

  fParticleGun->SetParticlePosition(pos);
  gParticleGun->SetParticlePosition(pos);

  fParticleGun->SetParticleMomentumDirection(mom);
  gParticleGun->SetParticleMomentumDirection(-1*mom);
  
  fParticleGun->SetParticleMomentum(511*keV);
  gParticleGun->SetParticleMomentum(511*keV);
  
    
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName="gamma";
  G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
  
  fParticleGun->SetParticleDefinition(particle);
  gParticleGun->SetParticleDefinition(particle);

  G4ThreeVector gammaMom = fParticleGun->GetParticleMomentumDirection();  
  G4ThreeVector pol1 = gammaMom.orthogonal().unit();  
  pol1 = pol1.rotate(2*CLHEP::pi * G4UniformRand(), gammaMom);     //randomizing the photon polarisation direction in the plane perp to gamma mom  
  fParticleGun->SetParticlePolarization(pol1);

  G4ThreeVector pol2 = pol1;
  pol2 = pol2.rotate(CLHEP::pi / 2.0, gammaMom);
  gParticleGun->SetParticlePolarization(pol2);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);    
  gParticleGun->GeneratePrimaryVertex(anEvent);
  
  // G4cout << "\n\nMomentum direction of the 1st incident photons: " << fParticleGun->GetParticleMomentumDirection () << G4endl;
  // G4cout << "polarisation vector (EF direction) of the 1st incident photons: " << fParticleGun->GetParticlePolarization() << "\n\n";

  // G4cout << "Momentum direction of the 2nd incident photons: " << gParticleGun->GetParticleMomentumDirection () << G4endl;
  // G4cout << "polarisation vector (EF direction) of the 2nd incident photons: " << gParticleGun->GetParticlePolarization() << "\n\n";

  // G4cout << "Angle betn pol1 and pol2: " << acos(pol1.dot(pol2))*(180/CLHEP::pi) << " degs\n\n";
}
