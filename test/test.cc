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
int main(int argc, char** argv)
{
  G4RunManager *runManager = new G4RunManager();
  runManager->SetUserInitialization(new MyDetectorConstruction());
  runManager->SetUserInitialization(new MyPhysicsList());
  runManager->SetUserInitialization(new MyActionInitialization());
  runManager->Initialize();

  G4UIExecutive *ui = 0;

  // G4EmParameters* params = G4EmParameters::Instance();
  // params->SetFluo(true);
  // params->SetAugerCascade(true);
  // params->SetAuger(false);
  // params->SetPixe(false);

  if(argc ==1)
    {
     ui = new G4UIExecutive(argc, argv);
    }
  
  G4VisManager *visManager = new G4VisExecutive();
  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  UImanager->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year"); // very high time threshold to allow all decays to happen
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle e- 0.01 mm");
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle e+ 0.01 mm");
  // UImanager->ApplyCommand("/run/setCutForAGivenParticle gamma 0.01 mm");
  UImanager->ApplyCommand("/run/initialize");          

  
  if(ui)
    {    
      UImanager->ApplyCommand("/control/execute vis.mac");
      //UImanager->ApplyCommand("/tracking/verbose 4");
      ui->SessionStart();
    }

  else
    {
      G4String command = "/control/execute ";
      //UImanager->ApplyCommand("/tracking/verbose 1");
      G4String fileName = argv[1];       
      UImanager->ApplyCommand(command+fileName);      
    }

  return 0;
}
