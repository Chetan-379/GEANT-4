#ifndef	GENERATOR_HH
#define	GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4IonTable.hh"
#include "G4Geantino.hh"
#include "G4StokesVector.hh"
#include "Randomize.hh"

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  MyPrimaryGenerator();
  ~MyPrimaryGenerator();

  virtual void GeneratePrimaries(G4Event*);

private:
  G4ParticleGun *fParticleGun;
};

#endif
