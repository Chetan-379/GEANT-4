#!/bin/bash

# Check if a ROOT file name is passed as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <ROOT file name>"
    exit 1
fi

root_file="$1"
material="$2"
cp test.root ../Analyze/src_root_files/"$1"
cd ../Analyze/
# Define the text file where you want to write the ROOT file name
txt_file="inputfile.txt"

# Clear the contents of the text file
> "$txt_file"

# Write the name of the ROOT file into the text file
echo "src_root_files/$1" > "$txt_file"

make

./analyzeLightBSM inputfile.txt out_root_files/"${root_file%.root}_out.root" $material
root -l "plot_all_hists.C(\"./\", \"${root_file%.root}_out.root\", \"$material\")" 
