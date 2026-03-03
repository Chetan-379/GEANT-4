#include <iostream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4EmParameters.hh"

#include "construction.hh"
#include "physics.hh"
#include "action.hh"

#include "Randomize.hh"
#include <chrono>
#include <unistd.h>   // getpid()

#include <string>

long gRandomSeed = 0;
std::string gOutputFileName = " ";

int main(int argc, char** argv)
{
   // -----------------------------
    // Argument handling
    // -----------------------------
    if (argc < 2) {
        G4cerr << "Usage: " << argv[0]
               << " run.mac [output.root]" << G4endl;
        //return 1;   //comment this if you want to run in the visualization
    }

    if (argc >= 3) {
      gOutputFileName = argv[2];
    }
    
  //-----------randomizing the run in terms of running time and system
  // High-resolution time
    auto now = std::chrono::high_resolution_clock::now();
    long timeSeed = now.time_since_epoch().count();

    // Combine with PID
    gRandomSeed = timeSeed ^ getpid();

    // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    // CLHEP::HepRandom::setTheSeed(gRandomSeed);
    
    G4cout << "\n###########################" << G4endl;
    G4cout << "Random seed = " << gRandomSeed << G4endl;
    G4cout << "OutFileName = " << gOutputFileName << G4endl;
    G4cout << "###########################" << G4endl;
    
  G4RunManager *runManager = new G4RunManager();
  runManager->SetUserInitialization(new MyDetectorConstruction());
  runManager->SetUserInitialization(new MyPhysicsList());
  runManager->SetUserInitialization(new MyActionInitialization());
  runManager->Initialize();

  G4UIExecutive *ui = 0;
  G4String gdmlFile = "";
  if(argc ==1)
    {
     ui = new G4UIExecutive(argc, argv);
    }

  // if(argc == 2)
  //   {
  //     ui = new G4UIExecutive(argc, argv);
  //     //gdmlFile = argv[1];
  //   }

  // gdmlFile = argv[1];
  
  G4VisManager *visManager = new G4VisExecutive();
  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  //UImanager->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year"); // very high time threshold to allow all decays to happen
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle e- 0.01 mm");
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle e+ 0.01 mm");
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle gamma 0.01 mm");
  // UImanager->ApplyCommand("/run/initialize");          

  
  if(ui)
    {    
      UImanager->ApplyCommand("/control/execute vis.mac");
      //UImanager->ApplyCommand("/tracking/verbose 1");
      ui->SessionStart();
    }

  else
    {
      G4String command = "/control/execute ";
      //UImanager->ApplyCommand("/tracking/verbose 1");
      G4String fileName = argv[1];       
      UImanager->ApplyCommand(command+fileName);

      G4RunManager::GetRunManager()->SetVerboseLevel(0);
    }

  return 0;
}
