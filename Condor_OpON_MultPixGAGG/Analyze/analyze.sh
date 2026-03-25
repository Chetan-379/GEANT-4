# make
# ./analyzeLightBSM inputfiles.txt loop_geom_theta_check.root GAGG
# root loop_geom_theta_check.root

make
./analyzeLightBSM inputfiles.txt loop_geom_phi_check_anyEinc.root GAGG
root loop_geom_phi_check_anyEinc.root
