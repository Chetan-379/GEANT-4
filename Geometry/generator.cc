#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
  //G4ThreeVector mom(1.5, 2.5, 0.);

  G4ThreeVector mom(1.5, 2.5, 1.);

  //G4ThreeVector mom(0, 1.5, 0);
  //G4StokesVector pol(G4ThreeVector(0.,0.,0.));
  
  //pol.SetPhoton();

  // G4cout << "\n\n\n\n\n=========== mom mag" << mom.mag() << G4endl;
  G4cout << " mom unit mag" << mom.unit()<< G4endl;
  
  //G4ThreeVector mom(1, 0, 0.);

  //G4ThreeVector mom(0, 2.5,1.5);
  //G4ThreeVector mom(0.,1., 0.);
  //G4ThreeVector mom(1.55, 2.5, 0.);
  
  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(511*keV);

  
  G4ThreeVector pol = fParticleGun->GetParticleMomentumDirection().orthogonal().unit();
  G4cout << "\n\npolarsation direction: " << pol << G4endl;
  fParticleGun->SetParticlePolarization(pol);
  //fParticleGun->SetParticlePolarization(G4ThreeVector(1,0,0));
  //fParticleGun->SetParticlePolarization(G4ThreeVector(1,0,0));

}

MyPrimaryGenerator::~MyPrimaryGenerator()
{
  delete fParticleGun;
}

void MyPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{    
    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName="gamma";
    G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
    
    fParticleGun->SetParticleDefinition(particle);  
    fParticleGun->GeneratePrimaryVertex(anEvent);

    auto pol = fParticleGun->GetParticlePolarization();
    G4cout << "\n\nMomentum direction of the incident photons: " << fParticleGun->GetParticleMomentumDirection () << G4endl;
    G4cout << "stokes vector of the incident photons: " << pol << "\n\n";
}
