- In this version, a machinery to run the Geant4 simulation on LXPLUS is added.
- Following is a typical sequence to run the simulation
  ```md
  git clone --branch G4_condor --depth 1 https://github.com/Chetan-379/GEANT-4.git
  ```
  
  ```md
  cd GEANT-4/Condor_OpON_MultPixGAGG/condor
  ```
   Do the relevant modifications in the paths and make the relevant directories in the relevant places.  
```md
python3 condor.py  
python3 submit.py  
```
