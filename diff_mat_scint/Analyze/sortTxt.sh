#!/bin/bash

# Check if directory argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

DIR="$1"
Mat="$2"

# Check if directory exists
if [ ! -d "$DIR" ]; then
    echo "Error: Directory '$DIR' not found!"
    exit 1
fi

# Loop through all .txt files in the directory
for f in "$DIR"/*.txt; do
    # Check if file exists (handles case of no .txt files)
    [ -e "$f" ] || { echo "No .txt files found in $DIR"; exit 0; }

    sort -n -k1,1 "$f" -o "$f"
done

mv $DIR/*.txt $DIR/$Mat/
