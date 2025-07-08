#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
  G4ThreeVector mom(1.5, 2.5, 0.);
  //G4ThreeVector mom(0.,1., 0.);
  //G4ThreeVector mom(1.55, 2.5, 0.);
  
  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(511*keV);
  //fParticleGun->SetParticlePolarization(G4ThreeVector(0,0,1));
  //fParticleGun->SetParticlePolarization(G4ThreeVector(1,0,0));
  fParticleGun->SetParticlePolarization(G4ThreeVector(0,0,1));

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
    G4cout << "\n\npolarisation of the photons: " << pol << G4endl;
}
