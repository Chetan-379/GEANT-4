In this version, mainly working machinery for calculating the angle between scattering plane and the polarisation plane (eta) is included.

- Here we have tried to reproduce the result from Klein Nishina formula for dependence of cross section on eta.
- We got the distribution of eta over 10000 events similar to the expected one (peak at +90 and -90 degrees).
- We took only first compton scattering of the event in the Analysis script based on the timing information of the hits.
- There was a confusion regarding `SetParticlePolarisation()`: whether it sets the direction of electric field or the stokes vector of the photon.
- We arrived at the conclusion that the polarisation vector there is the electric field direction. since simulation results are consistent only when the polarisation vector and photon momentum direction are perpendicular.
- As of now, we plan to move ahead with this interpretation, I will modify the things in future if any mistake is found.


=======To run the simulation:
- `mkdir build`
- `cmkae ..`
- `./test`