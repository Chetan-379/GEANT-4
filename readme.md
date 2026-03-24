In this tag, a scripts are added for the compton scattering angle reconstruction in a multipixelated CCal setup. In the analyzer a comparison of the reconstructed theta is done with the truth theta in the scatterer and the reco theta without the CCal Setup.

- right now the angle reconstruction is done using the true Edep info, not simulating the optical photons there.
- To simulate the optical photons uncomment the L146 and L147 in the construction.cc file.
- will have to write by your own the analyzer script to scattering angle reconstruction with optical photons (spread the chosen energy ranges according to the energy resolution).  

=================================  
To Run the simulation:  
`clone the git repository: git clone --branch CCal_OpOFF_angle_reco_multipixel --depth 1 https://github.com/Chetan-379/GEANT-4.git`

To simulate with CCal setup:  
`cp CCal_construction construction.cc`  
`mkdir build`  
`cd build`  
`cmake ..`  
`make`  
`./CCal run.mac`  
`cp test.root ../Analyze/src_root_files/OpOFF_CCal_100kEvts.root`  

To simulate without CCal setup:  
`cd ..`  
`cp noCCal_construction construction.cc`  
Repeat from step 2-6 of CCal case  
`cp test.root ../Analyze/src_root_files/OpOFF_noCCal_100kEvts.root`  

To do the analysis:  
`cd Analyze`  
`make`  
for CCal case:   
`./analyzeLightBSM CCal_inputfiles.txt OpOFF_CCal_100kEvts_out.root GAGG_LYSO`  
for no CCal case:   
`./analyzeLightBSM noCCal_inputfiles.txt OpOFF_noCCal_100kEvts_out.root GAGG`  

