#!/bin/bash

make
materials=("PbWO4" "BGO")

txt_file="inputfile.txt"
for id in "${materials[@]}"; do

    for f in src_root_files/*.root; do
	#> "$txt_file"
      [[ -e "$f" ]] || continue
      if [[ "$f" == *"$id"* ]]; then
	  #printf '%s\n' "$f" >> "$output_file"
	  echo "$f" > "$txt_file"

	  initfile="${f#src_root_files/}"
	  output_root_file="$out_root_files/${initfile%.root}_out.root"

	 echo $output_root_file

	 ./analyzeLightBSM inputfile.txt $output_root_file $id
	  	  
      fi
  done
done
