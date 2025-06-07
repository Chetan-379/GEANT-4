Steps to generate the Tree and Analyze the tree to calculate Resolution, Efficiency and make the plots for PbWO4 10cm thickness:

1. clone the git repository: git clone -b Efficiency_till_TIFR_Meet https://github.com/Chetan-379/GEANT-4.git
2. Go to the working directory: cd GEANT-4/diff_mat_scint
3. Generate the tree for 10 cm PbWO4: (for other thickness modify the line "G4double rad_dz = 10 * cm" in "PbWO4_construction construction.cc" file.)
   - cp PbWO4_construction construction.cc  (Use BGO instead of PbWO4 if you want to generate for BGO)
   - mkdir build && mkdir -p Analyze/src_root_files && mkdir -p Analyze/out_root_files && cd build (replace PbWO4 with BGO if want for BGO)
   - cmake ..
   - cp ../runJob.sh ./ && ./runJob.sh

4. Analyze the tree:
	- cd ../Analyze/ && ./make_hists.sh
	- ./make_hists.sh

5. Do the fitting and store the plots:
	- ./make_dir.sh
	- root -b ‘FitFunc.C(“./”)’ 
Plots can be found in the ‘plots’ directory

To make the plots of the physical variables of Scintillation photons use plot_all_hists.C
To overlay the data from different txt files use: Overlay_txt.C
To plot data from a single txt file use: plot_txt.C
