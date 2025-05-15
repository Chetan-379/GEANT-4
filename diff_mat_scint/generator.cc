#include "generator.hh"
MyPrimaryGenerator::MyPrimaryGenerator()
{
  fParticleGun = new G4ParticleGun(1);
  G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
  G4ThreeVector mom(0., 0., 1.);
  
  fParticleGun->SetParticlePosition(pos);
  fParticleGun->SetParticleMomentumDirection(mom);
  fParticleGun->SetParticleMomentum(0.5*MeV);
  //fParticleGun->SetParticleEnergy(0.0 * eV);
   
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
    
    // G4ThreeVector pos(0.*cm, 0.*cm, 0.*cm);
    // G4ThreeVector mom(0., 0., 0.);

    // if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {
    // fParticleGun->SetParticlePosition(pos);
    // fParticleGun->SetParticleMomentumDirection(mom);
    // fParticleGun->SetParticleMomentum(100.*GeV);
    fParticleGun->SetParticleDefinition(particle);

    // G4int Z = 9;
    // G4int A = 18;

    // G4int Z = 11;
    // G4int A = 22;

    // G4double charge = 0. * eplus;
    // G4double energy = 0 * keV;

    // G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(Z, A, energy);
    // fParticleGun->SetParticleDefinition(ion);
    // fParticleGun->SetParticleCharge(charge);
    // }
  

  //fParticleGun->SetParticleEnergy(energy);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //G4ParticleDefinition* na22 = G4ParticleTable::GetParticleTable()->FindParticle();
  // G4DecayTable* decayTable = ion->GetDecayTable();
  // G4cout << "\n\nDecayTable: " << decayTable << "\n\n" << G4endl; 
  // //if (decayTable) {
  //   G4int numChannels = decayTable->entries();
  //   G4cout << "Na-22 Decay Channels:" << G4endl;
  //   for (G4int i = 0; i < numChannels; ++i) {
  //     G4VDecayChannel* channel = decayTable->GetDecayChannel(i);
  //     G4cout << "Decay Channel " << i+1 << ": " << channel << G4endl;
  //     }
  //   //}

}
