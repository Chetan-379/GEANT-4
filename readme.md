In this version, further modifications are done for studying two back to back photons:
- hits are randomized in time.
- following changes in Event selection:
  - if a photon does multiple scattering in a single tile count them 1 only and take the information of first interaction.
  - post step point momentum at compton scattering point is added to calculate the angle between the scattering planes.
  - for nHits>=2 case, calculated the angle between the scattering planes in lower and higher visibility range of theta. 

=======To run the simulation:
- `mkdir build`
- `cmake ..`
- `./test`