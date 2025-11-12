In this version, have cleaned the scripts and further modifications are done for studying two back to back photons:
- hits are randomized in time.
- following changes in Event selection:
  - if a photon does multiple scattering in a single tile count them 1 only and take the information of first interaction.
  - post step point momentum at compton scattering point is added to calculate the angle between the scattering planes.
  - for nHits>=2 case, calculated the angle between the scattering planes in lower and higher visibility range of theta. 

=======To run the simulation:
- `mkdir build`
- `cmake ..`
- `./test`

========To run the Analyzer:

`cd Analyze`

`make`

`./analyzeLightBSM inputfile.txt outFile.root Plastic_BGO`

Note: inputfile.txt contains the name the root file containing the tree generated from the G4 Simulation. The root files mentioned there are merged together and a single output root file will be generated.