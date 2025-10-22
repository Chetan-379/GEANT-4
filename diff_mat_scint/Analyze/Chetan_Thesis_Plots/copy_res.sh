#!/bin/bash

# Set the base directory (the folder containing the "res" folder)
base_dir="./"

# Check if "res" directory exists in the base directory
if [[ ! -d "$base_dir/res" ]]; then
  echo "The 'res' folder does not exist in the base directory."
  exit 1
fi

# Loop through each subfolder in the base directory
#for subfolder in "$base_dir"/*/*/*/*/*; do
for subfolder in "$base_dir"*/*/*/*/*/*/*; do
  if [[ -d "$subfolder" && $(basename "$subfolder") != "res" ]]; then
    echo "Copying 'res' to $subfolder"
    cp -r "$base_dir/res" "$subfolder/"
  fi
done

echo "Done!"
