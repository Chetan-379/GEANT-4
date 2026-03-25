#!/bin/bash
#cd /afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana_FR

# eval `scram runtime -sh`

source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh
export X509_USER_PROXY=/afs/cern.ch/user/c/cagrawal/x509up_u165745

echo $PWD
cd /afs/cern.ch/user/c/cagrawal/public/run_G4/GEANT-4/OpON_MultPixGAGG
./build/$1 $2 $3

mv $3 /eos/home-c/cagrawal/G4_Sim_outputs/
mv $4 /eos/home-c/cagrawal/G4_Sim_outputs/logs/
mv $5 /eos/home-c/cagrawal/G4_Sim_outputs/logs/
mv $6 /eos/home-c/cagrawal/G4_Sim_outputs/logs/

