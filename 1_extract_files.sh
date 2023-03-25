#!/bin/bash

# Author: Tankred Ott
# Basic usage: 1_extract_files.sh [name_of_csv_file].csv [input_folder] [output_folder]

locus_file=$1
in_dir=$2
out_dir=$3

while IFS= read -r line; do
  echo "$line"
  ls -d ${in_dir}/* | grep -e "$line" | xargs -i mv {} $out_dir
done < "$locus_file"
