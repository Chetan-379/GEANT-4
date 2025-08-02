#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  //G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);

  //G4ThreeVector mom(1.5, 2.5, 0.);
  //G4ThreeVector mom(1.5, 2.5, 1.);

  //G4ThreeVector mom(0, 1.5, 0);
  //G4StokesVector pol(G4ThreeVector(-1., 0., 0.));
  //pol = G4ThreeVector(0.,1.,0.);
  
  // fParticleGun->SetParticlePosition(pos);
  // fParticleGun->SetParticleMomentumDirection(mom);
  // fParticleGun->SetParticleMomentum(511*keV);

  // G4ThreeVector gammaMom = fParticleGun->GetParticleMomentumDirection();
  // G4ThreeVector pol = gammaMom.orthogonal().unit();
  // pol = pol.rotate(2*CLHEP::pi * G4UniformRand(), gammaMom);     //randomizing the photon polarisation direction in the plane perp to gamma mom

  // G4cout << "\nrandom number: " << G4UniformRand() << G4endl;
  // G4cout << "\n\npolarsation direction: " << pol << G4endl;
  // fParticleGun->SetParticlePolarization(pol);
}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  fParticleGun = new G4ParticleGun(1);
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
  
  
  //G4ThreeVector mom(1.5, 2.5, 0.);
  G4ThreeVector mom(0, 1, 0.);
  G4ThreeVector ZAxis(0,0,1);
  G4ThreeVector XAxis(1,0,0);

  
  mom = mom.rotate(0.6044*((2*G4UniformRand()) - 1.0), XAxis);   //first randomizing in z
  mom = mom.rotate(2*CLHEP::pi * G4UniformRand(), ZAxis);        //then randomizing in xy plane 
  

  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(511*keV);
  
  G4ThreeVector gammaMom = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector pol = gammaMom.orthogonal().unit();
  
  pol = pol.rotate(2*CLHEP::pi * G4UniformRand(), gammaMom);     //randomizing the photon polarisation direction in the plane perp to gamma mom
  fParticleGun->SetParticlePolarization(pol);
  
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName="gamma";
  G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
  
  fParticleGun->SetParticleDefinition(particle);  
  fParticleGun->GeneratePrimaryVertex(anEvent);
  

  G4cout << "\n\nMomentum direction of the incident photons: " << fParticleGun->GetParticleMomentumDirection () << G4endl;
  G4cout << "polarisation vector (EF direction) of the incident photons: " << pol << "\n\n";
}
