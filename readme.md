In this version, Scripts are modified for studying two back to back photons:
- In the primary generator file another photon gun is added to shoot photons in the opposite direction to first one.
- All the hits in an event are analyzed after sorting them in the ascending order of time.
- In the Analyzer script a working machinery is added to calculate the scattering angle in a semirealistic way using the method mentioned in the JPET paper.

- IMPROVEMENTS NEEDED:
  - In this version of Analyzer, the 4 photon event yield necessary to calculate theta is too low. need to improve on that.
  - Randomizing the hits in time is needed.


=======To run the simulation:
- `mkdir build`
- `cmkae ..`
- `./test`