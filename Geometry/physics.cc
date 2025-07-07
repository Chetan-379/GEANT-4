#include "physics.hh"

MyPhysicsList::MyPhysicsList()
{
  //RegisterPhysics (new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4EmLivermorePolarizedPhysics());
  RegisterPhysics (new G4OpticalPhysics());
}

MyPhysicsList::~MyPhysicsList()
{}
