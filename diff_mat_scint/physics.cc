#include "physics.hh"

MyPhysicsList::MyPhysicsList()
{
  RegisterPhysics (new G4EmStandardPhysics_option4());
  //RegisterPhysics (new G4EmStandardPhysics());
  //RegisterPhysics (new G4EmLivermorePhysics());
  RegisterPhysics (new G4OpticalPhysics());
  // RegisterPhysics (new G4RadioactiveDecayPhysics());
  // RegisterPhysics (new G4DecayPhysics());
}

MyPhysicsList::~MyPhysicsList()
{}
