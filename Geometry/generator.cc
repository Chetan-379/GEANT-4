#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
  G4ThreeVector mom(1.5, 2.5, 0.);
  //G4ThreeVector mom(1.55, 2.5, 0.);
  
  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  //fParticleGun->SetParticleMomentum(0.5*MeV);
  fParticleGun->SetParticleMomentum(511*keV);   
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
}
