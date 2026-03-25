#!/bin/bash
make
Energies=(511 100 150 300 450 600 800 1000)

for energy in "${Energies[@]}"; do
    file="test.mac"
    echo "/gun/momentumAmp $energy keV" >> "$file"
    echo "/run/beamOn 10000" >> "$file"

    ./test test.mac
    cp test.root ../Analyze/src_root_files/PbWO4_10cm_ScintON_${energy}keV_Birk0_NoAbs_NoCeren_DerEmis_spline.root
    echo Energy is: $energy keV
    echo test.root shifted to ../Analyze/src_root_files/PbWO4_10cm_ScintON_${energy}keV_Birk0_NoAbs_NoCeren_DerEmis_spline.root
     : > "$file"

done
