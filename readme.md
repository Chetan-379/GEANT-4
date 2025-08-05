In this version, most of the changes are made in the Analyzer script and the Primary generator script.

-Definition of nHit variable (to calculate no. of Hits in the event) is updated to include only the events where 1st hit is from Compton Scattering.
-Added in the Analyzer script to analyze all the events into different categories based on location of hits. e.g 1st hit in plastic 2nd in BGO etc.

-In the primary generator file, randomized the incident photon momentum direction in the XY plane and also along the axis of the cylinder.
-Also randomized the direction of polarisation of incident photon in the plane perpendicular to the momentum direction.

Note: The 2d distribution of Analzer theta Vs step theta won't behave as expected. To see the expected behaviour change the 'nHitComp>0' to 'nHitComp>1' in the Analyzer script.

=======To run the simulation:
- `mkdir build`
- `cmkae ..`
- `./test`